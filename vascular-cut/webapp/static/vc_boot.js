// Small things the page needs that trame has no widget for.
//
// Loaded as a module script beside vc_upload.js: server.enable_module({"scripts": ...}).
// Templates reach these through utils.get("vcCopy") and so on; a template that calls
// one during render guards against it not being loaded yet.

// Copy text to the clipboard. The submission ID is the only way to have uploaded
// files found and removed, so copying it has to work: the clipboard API first, then
// the old select-and-copy for a page that is not a secure context.
window.vcCopy = function (text) {
  text = String(text || "");
  if (!text) return;
  function fallback() {
    try {
      var ta = document.createElement("textarea");
      ta.value = text;
      ta.setAttribute("readonly", "");
      ta.style.position = "fixed";
      ta.style.opacity = "0";
      document.body.appendChild(ta);
      ta.select();
      document.execCommand("copy");
      document.body.removeChild(ta);
    } catch (e) { /* nothing more to try; the ID is still on screen */ }
  }
  try {
    if (navigator.clipboard && window.isSecureContext) {
      navigator.clipboard.writeText(text).catch(fallback);
      return;
    }
  } catch (e) { /* fall through */ }
  fallback();
};

// Scroll something in the drawer into view, once whatever opened it has rendered
// (a panel expanding, a card growing its details). Two frames, then a short wait for
// the panel's expand transition.
window.vcReveal = function (selector) {
  function go() {
    var el = document.querySelector(selector);
    if (el) el.scrollIntoView({ block: "start", behavior: "smooth" });
  }
  requestAnimationFrame(function () {
    requestAnimationFrame(function () { setTimeout(go, 260); });
  });
};

// The same, for one opening's card by its label.
window.vcRevealOp = function (label) {
  if (label == null) return;
  var esc = window.CSS && CSS.escape ? CSS.escape(String(label)) : String(label);
  window.vcReveal('.op-card[data-label="' + esc + '"]');
};

// Why an opening will not cut, in the server's own words: the status line opens with
// the opening's name when an opening is marked bad, and its first sentence is the
// reason. Nothing is reworded; "" when the status is about something else.
window.vcWhy = function (status, name) {
  status = String(status || "");
  name = String(name || "");
  if (!name || status.indexOf(name) !== 0) return "";
  var rest = status.slice(name.length).replace(/^[:,\s]+/, "");
  var first = rest.split(". ")[0].replace(/\.$/, "");
  if (!first) return "";
  return first.charAt(0).toUpperCase() + first.slice(1) + ".";
};

// How many rim slivers the status line reports, for the Remove rim slivers row.
window.vcSlivers = function (status) {
  var m = /(\d+) rim sliver/.exec(String(status || ""));
  return m ? m[1] + " on the rims now" : "";
};

// A path field shows the end of its path (the file name) whenever it is not being
// edited: the start is a folder every case shares.
(function () {
  var SEL = '[data-test="path-mesh"] input, [data-test="path-centerline"] input, '
          + '[data-test="sm-path"] input';
  function tail(el) { if (el && document.activeElement !== el) el.scrollLeft = el.scrollWidth; }
  document.addEventListener("focusout", function (e) {
    if (e.target && e.target.matches && e.target.matches(SEL)) {
      setTimeout(function () { tail(e.target); }, 0);
    }
  });
  setInterval(function () {
    document.querySelectorAll(SEL).forEach(tail);
  }, 1500);
})();

// The file name from a path, for the phone, where a path field shows only its
// folder.
window.vcBase = function (path) {
  return String(path || "").split(/[\\/]/).pop();
};

// A phone opens on the model, not on a drawer covering it. Below 960 px the drawer
// is an overlay; start it closed. Waits for the server's first state, which would
// otherwise arrive afterwards and open it again.
(function () {
  var tries = 0;
  function start() {
    var t = window.trame;
    var ready = t && t.state && t.state.get && t.state.get("mode") !== undefined;
    if (!ready) {
      if (tries++ < 100) setTimeout(start, 100);
      return;
    }
    if (window.innerWidth < 960) {
      t.state.set("main_drawer", false);
    }
  }
  start();
})();
