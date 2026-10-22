// Uploads and downloads over plain HTTP, beside trame's websocket rather than through it.
//
// Everything trame does travels one websocket, in order, and its client holds back every
// state update while a message is in flight. So a file sent that way freezes the whole
// page -- the 3-D view, every control, every message -- for as long as it takes to move:
// measured at about 10 s for 20 MB at 2 MiB/s, with no renders and no response at all.
// An HTTP request on the same server runs alongside the socket, so the page stays live,
// and the browser can report progress as the bytes go. The server side is the pair of
// routes cut_app.py adds through on_server_bind.
//
// Loaded as a module script: server.enable_module({"scripts": [...]}). Templates reach
// these through utils.get("vcUpload") / utils.get("vcDownload").

// POST one File to url. Resolves {status, body}; never rejects, so a template can chain
// on it without an unhandled rejection. Progress (0-100) goes into trame state.
window.vcUpload = function (file, url, token, progressKey) {
  return new Promise(function (resolve) {
    if (!file) { resolve(null); return; }
    var xhr = new XMLHttpRequest();
    xhr.open("POST", url);
    xhr.setRequestHeader("X-VC-Token", token);
    // The name travels as a header, URL-encoded. It is metadata only: the server never
    // uses it as a path, and stores the file under a generated id.
    xhr.setRequestHeader("X-VC-Name", encodeURIComponent(file.name));
    xhr.setRequestHeader("Content-Type", "application/octet-stream");
    var last = 0;
    if (progressKey) window.trame.state.set(progressKey, 0);
    xhr.upload.onprogress = function (e) {
      var now = Date.now();
      if (progressKey && e.lengthComputable && now - last > 200) {
        last = now;
        window.trame.state.set(progressKey, Math.round(100 * e.loaded / e.total));
      }
    };
    xhr.onload = function () {
      // Only this app's JSON carries a message worth showing. Anything else -- a
      // proxy's HTML error page, say -- becomes an empty body, and the page words it.
      var body = {};
      try { body = JSON.parse(xhr.responseText); } catch (err) { body = {}; }
      if (!body || typeof body !== "object") body = {};
      if (typeof body.error !== "string") body.error = "";
      if (progressKey) window.trame.state.set(progressKey, -1);
      resolve({ status: xhr.status, body: body });
    };
    xhr.onerror = function () {
      if (progressKey) window.trame.state.set(progressKey, -1);
      resolve({ status: 0, body: { error: "The connection dropped during the upload." } });
    };
    xhr.send(file);
  });
};

// Follow a one-time download link the server has just minted. An <a download> click
// rather than window.location, so the page is never navigated away from.
window.vcDownload = function (href) {
  if (!href) return;
  var a = document.createElement("a");
  a.href = href;
  a.download = "";
  a.rel = "noopener";
  document.body.appendChild(a);
  a.click();
  a.remove();
};
