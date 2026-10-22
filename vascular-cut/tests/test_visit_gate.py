"""The visit gate admits one browser, refuses the rest before they get any state, and
reports the admitted one leaving -- and nothing else.

No server, no browser: webapp/visit_gate.py is driven directly, with stand-ins for
the two wslink objects it touches. The browser-level proof, against the real app, is
tests/web/check_privacy.py; this is the fast one that pins each rule on its own.

WHAT IT PINS, each a way the privacy fix could quietly stop working:
  1. A claim hands out a ticket; while that ticket or an admitted visitor stands, every
     other claim gets the in-use page (409) and no ticket.
  2. Only the current ticket passes the hello -- not trame's public default secret,
     not a wrong one, not the same ticket twice, not an expired one.
  3. A refused socket closing does not count as the visitor leaving (it would wipe
     the person who IS connected); the admitted one closing does, exactly once.
  4. Every admit and every leave changes gen, which is what work in flight checks.
  5. If the app cannot make the page clean on admit, the visitor is refused.
  6. A reload -- the old socket closing while the new page is already claiming --
     lets the same person back in rather than refusing them as someone else.
  7. A tab that sends nothing past the idle limit has its socket closed; one that
     keeps sending does not.
  8. Wiring: the claim route is added, and wslink's application.exit RPCs, which let
     any admitted client stop the process, are removed.
  9. A page may write only the state keys the app owns. trame accepts any key from an
     admitted page -- its own trame__template_main, trame__scripts, trame__vue_use
     among them, which the NEXT visitor's page is built from and runs -- so without
     this one visitor could plant script in the page the next one uploads through.
     Lifecycle calls are limited to the ones trame's page sends, and an RPC that
     raises no longer ships the server's traceback, with its paths, to the page.
"""
import asyncio
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "webapp"))
from visit_gate import VisitGate                                   # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        fails.append(name)


class Calls:
    def __init__(self):
        self.admit = self.leave = 0

    def on_admit(self):
        self.admit += 1

    def on_leave(self):
        self.leave += 1


def make(**kw):
    c = Calls()
    g = VisitGate(busy_html="<p>in use</p>", on_admit=c.on_admit, on_leave=c.on_leave,
                  settle=kw.pop("settle", 0.2), log=lambda *_: None, **kw)
    return g, c


def body(resp):
    return json.loads(resp.text)


async def main():
    print("=" * 78)
    print("1-4. claim, admit, refuse, leave")
    print("=" * 78)
    g, c = make()
    r = await g.claim()
    t1 = body(r).get("secret")
    check("a claim on an empty instance gets 200 and a ticket",
          r.status == 200 and isinstance(t1, str) and len(t1) > 20, str(r.status))
    check("  ...with the session URL the page will open",
          body(r).get("sessionURL") == "ws://USE_HOST/ws")
    check("  ...and is never cached", r.headers.get("Cache-Control") == "no-store")
    r2 = await g.claim()
    check("a second claim while that ticket stands gets the in-use page, no ticket",
          r2.status == 409 and body(r2) == {"error": "<p>in use</p>"}, str(r2.status))
    check("trame's public default secret is refused",
          not await g.validate("wslink-secret", "raw"))
    check("a wrong token is refused", not await g.validate("x" * 32, "wrong"))
    check("a non-string token is refused", not await g.validate(None, "none"))
    check("nothing admitted yet", g.holder is None and c.admit == 0)
    gen0 = g.gen
    check("the ticket admits its page", await g.validate(t1, "A") and g.holder == "A")
    check("  ...the app was told once, and gen moved", c.admit == 1 and g.gen != gen0)
    check("the same ticket cannot be used twice", not await g.validate(t1, "A2"))
    t0 = time.monotonic()
    r3 = await g.claim()
    waited = time.monotonic() - t0
    check("a claim while someone is admitted gets the in-use page",
          r3.status == 409 and "secret" not in body(r3), str(r3.status))
    check("  ...after waiting the settle time for a reload to clear",
          0.15 <= waited < 1.0, f"{waited:.2f} s")
    gen1 = g.gen
    g.closed("some-refused-socket")
    g.closed(None)
    check("a refused socket closing is not the visitor leaving",
          g.holder == "A" and c.leave == 0 and g.gen == gen1)
    g.closed("A")
    check("the admitted socket closing is the visitor leaving: told once, gen moved",
          g.holder is None and c.leave == 1 and g.gen != gen1)
    g.closed("A")
    check("  ...and only once", c.leave == 1)
    r4 = await g.claim()
    t4 = body(r4).get("secret")
    check("after they leave, the next person gets a new ticket",
          r4.status == 200 and t4 and t4 != t1)
    check("  ...and is admitted", await g.validate(t4, "B") and c.admit == 2)

    print("\n" + "=" * 78)
    print("2. tickets expire")
    print("=" * 78)
    g, c = make(ticket_ttl=0.1)
    t = body(await g.claim())["secret"]
    await asyncio.sleep(0.15)
    check("an expired ticket is refused", not await g.validate(t, "late"))
    r = await g.claim()
    check("an expired, unused ticket no longer blocks the next claim", r.status == 200)

    print("\n" + "=" * 78)
    print("5. fail closed")
    print("=" * 78)

    def boom():
        raise RuntimeError("scene would not clear")
    g = VisitGate(busy_html="x", on_admit=boom, on_leave=lambda: None, settle=0.1,
                  log=lambda *_: None)
    t = body(await g.claim())["secret"]
    ok = await g.validate(t, "C")
    check("when the page cannot be made clean, the visitor is refused",
          ok is False and g.holder is None)
    check("  ...and the instance is free again, not stuck",
          (await g.claim()).status == 200)

    print("\n" + "=" * 78)
    print("6. reload: the old socket closes while the new page is claiming")
    print("=" * 78)
    g, c = make(settle=2.0)
    await g.validate(body(await g.claim())["secret"], "M-old")
    task = asyncio.ensure_future(g.claim())
    await asyncio.sleep(0.1)
    check("the new page's claim is waiting, not refused", not task.done())
    g.closed("M-old")
    r = await task
    check("once the old socket closes, the same person gets back in",
          r.status == 200 and await g.validate(body(r)["secret"], "M-new"))
    check("  ...and their old visit was wiped in between", c.leave == 1 and c.admit == 2)

    print("\n" + "=" * 78)
    print("7-8. idle limit, and the wiring")
    print("=" * 78)

    class FakeWS:
        def __init__(self):
            self.closed_with = None

        async def close(self, code=None, message=b""):
            self.closed_with = code

    calls = {"update": [], "life": [], "err": []}

    def real_update(obj, changes):
        calls["update"].append(changes)
        return "updated"

    def real_life(obj, name):
        calls["life"].append(name)

    class FakeHandler:
        def __init__(self):
            self.functionMap = {"application.exit": 1, "application.exit.later": 2,
                                "trame.state.get": 3,
                                "trame.state.update": ("core", real_update),
                                "trame.lifecycle.update": ("core", real_life)}
            self.connections = {}
            self.seen = []

        async def onCompleteMessage(self, rpc, client_id):
            self.seen.append(client_id)

        async def sendWrappedError(self, rpcid, code, message, data=None,
                                   client_id=None):
            calls["err"].append(data)

    class FakeRouter:
        def __init__(self):
            self.posts = []

        def add_post(self, path, fn):
            self.posts.append(path)

    class FakeServer:
        def __init__(self):
            self.app = type("A", (), {"router": FakeRouter()})()
            self._ws_handlers = [FakeHandler()]

    srv = FakeServer()
    h = srv._ws_handlers[0]
    g, c = make(idle_limit=1.5, settle=0.1,
                may_write=lambda k: not k.startswith("trame__") and k in {"opacity"})
    g.attach(srv)
    check("the claim route is added", "/paraview/" in srv.app.router.posts)
    check("application.exit and application.exit.later are gone",
          "application.exit" not in h.functionMap
          and "application.exit.later" not in h.functionMap
          and "trame.state.get" in h.functionMap)

    obj, upd = h.functionMap["trame.state.update"]
    r = upd(obj, [{"key": "opacity", "value": 10},
                  {"key": "trame__template_main", "value": "<script>"},
                  {"key": "trame__vue_use", "value": ["window.eval"]},
                  {"key": "invented", "value": 1}])
    check("a page's write to a key the app owns goes through",
          calls["update"] == [[{"key": "opacity", "value": 10}]] and r == "updated")
    check("  ...its writes to trame's own keys and to invented keys do not",
          all(ch["key"] == "opacity" for batch in calls["update"] for ch in batch))
    upd(obj, [{"key": "trame__scripts", "value": ["https://x.example/a.js"]}])
    check("a batch of only refused writes never reaches trame", len(calls["update"]) == 1)
    upd(obj, "not a list")
    check("  ...nor does a malformed one", len(calls["update"]) == 1)
    obj, life = h.functionMap["trame.lifecycle.update"]
    life(obj, "client_exited")
    life(obj, "server_bind")
    life(obj, "server_ready")
    check("lifecycle calls: trame's page's own pass, anything else does not",
          calls["life"] == ["client_exited"], f"{calls['life']}")
    await h.sendWrappedError(7, -32001, "Exception raised",
                             {"method": "trame.trigger", "exception": "KeyError('x')",
                              "trace": "File C:\\Users\\someone\\app.py, line 1"})
    check("an RPC that raised tells the page which call, not the traceback",
          calls["err"] == [{"method": "trame.trigger"}], f"{calls['err']}")
    await h.sendWrappedError(8, -32000, "Authentication failed")
    check("  ...and other errors pass unchanged", calls["err"][-1] is None)
    await g.validate(body(await g.claim())["secret"], "I")
    ws = FakeWS()
    h.connections["I"] = ws
    for _ in range(4):                       # keep sending for ~2.4 s
        await asyncio.sleep(0.6)
        await h.onCompleteMessage({"id": 1}, "I")
    check("messages still reach wslink through the activity counter", h.seen[-1] == "I")
    check("a tab that keeps sending is not closed", ws.closed_with is None)
    await asyncio.sleep(3.0)                 # now silent past the 1.5 s limit
    check("a tab silent past the idle limit has its socket closed",
          ws.closed_with == 4001, str(ws.closed_with))

    print("\n" + "=" * 78)
    print("the wslink face")
    print("=" * 78)
    g, c = make()
    p = g.protocol()
    t = body(await g.claim())["secret"]
    check("validateToken goes through the gate",
          await p.validateToken(t, "P") and g.holder == "P")
    p.onClose("P")
    check("onClose goes through the gate", g.holder is None and c.leave == 1)
    check("it exports no RPCs of its own",
          not [k for k in dir(p) if getattr(getattr(p, k, None), "_wslinkuris", None)])


asyncio.run(main())
print()
if fails:
    print(f"{len(fails)} FAILURE(S):")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("ALL PASS")
