"""One visitor per process, and nothing of theirs left for the next one.

WHY THIS EXISTS. trame keeps ONE state per process and pushes every change to every
connected browser, and a browser that connects is handed the whole current state
before anything else. So two people on one instance share a case, a cut, a camera and
a submission id, and someone arriving after another has left -- within the idle
shutdown, 300 s -- is handed the last person's case. Uploads are kept and may be
patient-derived, so that is a privacy failure, not a quirk.

WHAT THIS DOES. It admits one browser at a time, refuses everyone else BEFORE they
receive a byte of state, and tells the app the moment the admitted browser's socket
closes, so the app can wipe everything before the next person is let in.

HOW, in the installed trame 3.13.2 / trame-client 3.13.2 / wslink 2.5.7 -- pinned in
requirements.txt, because two of the three steps lean on their internals:

  1. CLAIM. Before opening its websocket, the trame page POSTs to paraview/ asking
     where to connect (its launcher protocol). This answers that request: 409 with a
     static "in use" page when someone holds the instance -- the page shows exactly
     that text and never opens a socket -- or 200 with a one-time ticket as the
     socket's secret.
  2. ADMIT. wslink awaits validateToken on every link protocol during the hello, and
     an unauthenticated socket can call nothing, fetch no state and receives no
     publishes. Only the current ticket passes, once, and only while nobody holds the
     instance. The check and the claim have no await between them, so two hellos
     cannot both get in. A client that skips step 1 and presents trame's public
     default secret is refused here: once a link protocol defines validateToken,
     wslink no longer compares the secret itself.
  3. LEAVE. wslink calls onClose from the socket handler's `finally`, for every way a
     socket ends: tab closed, reload, crash, a dead network once the heartbeat gives
     up. trame's own client_exited is sent only from beforeunload and cannot be relied
     on. onClose also fires for refused sockets, so only the holder's close counts.

What a page may do once admitted is narrowed too, because what one visitor leaves in
the process is what the next is served: it may write only the state keys the app owns
(never trame's own, which carry the page's template and script lists), call only the
lifecycle hooks trame's page sends, not stop the process, and not receive server
tracebacks.

Every admit and every leave bumps `gen`. Work started for a visitor (an upload still
streaming, a session still building in a thread) captures gen and drops its result if
it changed, so it cannot land on the next person's page.

The app's side -- what "wipe" means -- is not here. This module only decides who is in.
"""
import asyncio
import secrets
import time

from aiohttp import web
from wslink.websocket import LinkProtocol

NO_STORE = {"Cache-Control": "no-store"}


LIFECYCLE_OK = ("client_connected", "client_exited", "client_unmounted")


class VisitGate:
    def __init__(self, *, busy_html, on_admit, on_leave, may_write=None,
                 session_url="ws://USE_HOST/ws", ticket_ttl=10.0, settle=3.0,
                 idle_limit=0.0, clock=time.monotonic, log=print):
        # Rendered by the refused page with innerHTML, so it must be a constant: never
        # put anything a visitor supplied into it.
        self.busy_html = busy_html
        self.on_admit, self.on_leave = on_admit, on_leave
        # Which state keys a browser may write (a predicate on the key). trame accepts
        # ANY key from an admitted client -- including its own trame__template_main,
        # trame__scripts and trame__vue_use, which the NEXT visitor's page builds
        # itself from and runs. Unfiltered, one visitor could plant script in the
        # page the next one uploads their files through. None allows everything.
        self.may_write = may_write
        # The client lowercases this whole string and swaps USE_HOST for the page's
        # host, and ws: for wss: under https. A deployment under a path prefix needs
        # the prefix here, in lower case.
        self.session_url = session_url
        # How long a claim reserves the instance for the page that made it. The page
        # opens its socket within a second of claiming; this is generous for a slow
        # one, and short so that a claim nobody redeems locks nobody out for long.
        self.ticket_ttl = float(ticket_ttl)
        # A reload closes the old socket and claims again within milliseconds, in
        # either order. Waiting briefly for the holder to go lets the same person back
        # in rather than refusing them as "someone else".
        self.settle = float(settle)
        # An open but untouched tab keeps its socket alive indefinitely (browsers
        # answer pings), which would hold the instance forever. 0 turns this off.
        self.idle_limit = float(idle_limit)
        self.clock, self.log = clock, log

        self.holder = None            # the admitted socket's client_id
        self.ticket = None
        self.ticket_until = 0.0
        self.gen = 0
        self.last_active = 0.0
        self._handler = None          # wslink's socket handler, once attached
        self._idle_task = None

    # ── wiring ────────────────────────────────────────────────────────────────
    def protocol(self):
        """The link protocol to register with trame's add_protocol_to_configure."""
        return _GateProtocol(self)

    def attach(self, wslink_server):
        """At on_server_bind, while aiohttp's router is still open for additions."""
        app = wslink_server.app
        for path in ("/paraview/", "/paraview"):
            app.router.add_post(path, self.claim)
        handler = (getattr(wslink_server, "_ws_handlers", None) or [None])[0]
        self._handler = handler
        if handler is None:
            self.log("[visit] no wslink handler found: idle limit and exit lock off")
            self.idle_limit = 0.0
            return
        # Any admitted client could otherwise stop the process with these, taking the
        # instance down for whoever comes next.
        fm = getattr(handler, "functionMap", {})
        fm.pop("application.exit", None)
        fm.pop("application.exit.later", None)
        self._filter_rpcs(fm)
        # An exception in any RPC goes back to the page with its full traceback --
        # absolute paths, the account name. The server log keeps it; the page gets
        # only which call failed.
        send_err = getattr(handler, "sendWrappedError", None)
        if send_err is not None:
            async def quiet(rpcid, code, message, data=None, client_id=None):
                if isinstance(data, dict) and "trace" in data:
                    data = {"method": data.get("method")}
                return await send_err(rpcid, code, message, data, client_id=client_id)
            handler.sendWrappedError = quiet
        # Activity, for the idle limit: every complete message the holder sends.
        orig = getattr(handler, "onCompleteMessage", None)
        if orig is None:
            self.log("[visit] wslink has no onCompleteMessage: idle limit off")
            self.idle_limit = 0.0
            return

        async def counted(rpc, client_id):
            self.touch(client_id)
            return await orig(rpc, client_id)
        handler.onCompleteMessage = counted

    def _filter_rpcs(self, fm):
        """Put a filter in front of the two trame RPCs a page could abuse."""
        entry = fm.get("trame.state.update")
        if entry is not None and self.may_write is not None:
            obj, update = entry

            def update_state(o, changes, *a, **kw):
                changes = changes if isinstance(changes, list) else []
                kept = [c for c in changes if isinstance(c, dict)
                        and isinstance(c.get("key"), str) and self.may_write(c["key"])]
                if len(kept) != len(changes):
                    names = sorted({repr(str(c.get("key"))[:40]) for c in changes
                                    if isinstance(c, dict) and c not in kept})
                    self.log(f"[visit] refused a page write to {len(changes) - len(kept)}"
                             f" state key(s) it does not own: {', '.join(names[:5])}")
                return update(o, kept, *a, **kw) if kept else None
            fm["trame.state.update"] = (obj, update_state)
        # trame.lifecycle.update calls ANY ctrl.on_<name> the page names, with no
        # arguments -- on_server_bind among them. The page only ever sends these.
        entry = fm.get("trame.lifecycle.update")
        if entry is not None:
            obj, life = entry

            def life_cycle_update(o, name, *a, **kw):
                if name not in LIFECYCLE_OK:
                    self.log(f"[visit] refused lifecycle call {str(name)[:40]!r}")
                    return None
                return life(o, name, *a, **kw)
            fm["trame.lifecycle.update"] = (obj, life_cycle_update)

    # ── the three steps ───────────────────────────────────────────────────────
    async def claim(self, request=None):
        """POST paraview/: a ticket for this page, or the in-use page."""
        await self._settle()
        now = self.clock()
        if self.holder is not None or (self.ticket and now < self.ticket_until):
            return web.json_response({"error": self.busy_html}, status=409,
                                     headers=NO_STORE)
        self.ticket = secrets.token_urlsafe(24)
        self.ticket_until = now + self.ticket_ttl
        return web.json_response({"sessionURL": self.session_url,
                                  "secret": self.ticket}, headers=NO_STORE)

    async def validate(self, token, client_id):
        """The hello. True admits this socket; anything else is refused."""
        t = self.ticket
        if not (t and isinstance(token, str) and self.clock() < self.ticket_until
                and secrets.compare_digest(token.encode("utf-8"), t.encode("utf-8"))):
            return False
        if self.holder is not None:
            return False
        # No await from here to the end: the check above and this claim are one step.
        self.holder, self.ticket = client_id, None
        self.gen += 1
        self.last_active = self.clock()
        try:
            self.on_admit()
        except Exception as e:                                   # noqa: BLE001
            # Fail closed. If the page could not be made clean, nobody sees it.
            self.log(f"[visit] admit failed, refusing: {type(e).__name__}: {e}")
            self.holder = None
            self.gen += 1
            return False
        self.log(f"[visit] admitted (gen {self.gen})")
        if self.idle_limit and self._idle_task is None:
            self._idle_task = asyncio.ensure_future(self._watch_idle())
        return True

    def closed(self, client_id):
        """A socket ended. Only the holder's ending is a visitor leaving."""
        if client_id is None or client_id != self.holder:
            return
        self.holder = None
        self.gen += 1
        try:
            self.on_leave()
            self.log(f"[visit] ended, wiped (gen {self.gen})")
        except Exception as e:                                   # noqa: BLE001
            # The next admit wipes again before letting anyone in, and refuses them if
            # that fails too -- so a failure here is logged, not fatal.
            self.log(f"[visit] wipe on leave failed: {type(e).__name__}: {e}")

    # ── support ───────────────────────────────────────────────────────────────
    def touch(self, client_id):
        if client_id is not None and client_id == self.holder:
            self.last_active = self.clock()

    async def _settle(self):
        end = self.clock() + self.settle
        while self.holder is not None and self.clock() < end:
            await asyncio.sleep(0.05)

    async def _watch_idle(self):
        """Close the holder's socket once it has sent nothing for idle_limit seconds.

        Closing it is the whole eviction: wslink then runs onClose, which wipes. The
        page shows trame's "Connection closed", and the next visitor can get in.
        """
        while self.idle_limit:
            await asyncio.sleep(max(1.0, min(30.0, self.idle_limit / 4)))
            if self.holder is None or self._handler is None:
                continue
            if self.clock() - self.last_active < self.idle_limit:
                continue
            ws = getattr(self._handler, "connections", {}).get(self.holder)
            if ws is None:
                continue
            self.log(f"[visit] idle for {self.idle_limit:.0f} s, closing")
            try:
                await ws.close(code=4001, message=b"idle")
            except Exception as e:                               # noqa: BLE001
                self.log(f"[visit] could not close the idle socket: {e}")


class _GateProtocol(LinkProtocol):
    """The wslink face of the gate. Exports no RPCs."""

    def __init__(self, gate):
        super().__init__()
        self._gate = gate

    async def validateToken(self, token, client_id):
        return await self._gate.validate(token, client_id)

    def onClose(self, client_id):
        self._gate.closed(client_id)
