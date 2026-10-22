"""Browser regression for cut sliders and the two linked WebGL cameras.

Run with: .conda-env/bin/python -m unittest discover -s code/tests -v
Requires Playwright and its Chromium browser. API fixtures avoid loading or
modifying patient data; the page, Plotly, mouse drags, and slider handlers are real.
"""

import json
from pathlib import Path
import sys
import unittest
from urllib.parse import urlparse

from playwright.sync_api import TimeoutError as BrowserTimeout, sync_playwright
from plotly.offline import get_plotlyjs

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from review_cardiac_split_web import PAGE


def mesh(scale, shift=0):
    vertices = [(scale + shift, 0, 0), (-scale + shift, 0, 0),
                (shift, scale, 0), (shift, -scale, 0),
                (shift, 0, scale), (shift, 0, -scale)]
    faces = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
             (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    return {**dict(zip(('x', 'y', 'z'), map(list, zip(*vertices)))),
            **dict(zip(('i', 'j', 'k'), map(list, zip(*faces)))),
            'n_faces': len(faces)}


EMPTY = dict(x=[], y=[], z=[], i=[], j=[], k=[], n_faces=0)
QUALITY = dict(components=1, simple_boundary_loops=True, faces=8,
               boundary_loops=1, euler=1)
SNAPSHOT = """() => Object.fromEntries(['inner','outer'].map(id => {
    const gd = byId(id), scene = gd._fullLayout.scene;
    return [id, {camera: scene._scene.getCamera(),
        ranges: ['xaxis','yaxis','zaxis'].map(key => scene[key].range.slice()),
        aspect: {...scene.aspectratio},
        plane: gd.data.find(trace => trace.color === '#e63946').x.slice(),
        mesh: gd.data[0].x.slice()}];
}))"""


class ReviewCameraTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.playwright = sync_playwright().start()
        cls.browser = cls.playwright.chromium.launch(headless=True, args=[
            '--enable-webgl', '--use-gl=angle', '--use-angle=swiftshader',
            '--enable-unsafe-swiftshader',
        ])
        cls.plotly_js = get_plotlyjs()

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()

    def route(self, route):
        path = urlparse(route.request.url).path
        if path == '/':
            return route.fulfill(body=PAGE, content_type='text/html')
        if path == '/plotly.min.js':
            return route.fulfill(body=self.plotly_js,
                                 content_type='application/javascript')
        if path == '/api/cases':
            data = dict(cases=['camera-regression'], source_mode='raw_dataset')
        elif path == '/api/mesh':
            data = dict(case='camera-regression', inner=mesh(8), source=mesh(10),
                        bounds=[-10, 10, -10, 10, -10, 10],
                        plane=dict(origin=[3, 0, 0], normal=[1, 0, 0],
                                   axis_center=[0, 0, 0], axis_range=[-10, 10],
                                   method='fixture'), apply_cut_to_inner=True)
        elif path == '/api/preview':
            request = route.request.post_data_json
            shift = (request['origin'][0] * .01 + request['normal'][1] * .03
                     + request['normal'][2] * .02)
            data = dict(inner_kept=mesh(8, shift), inner_removed=EMPTY,
                        kept=mesh(10, shift), removed=EMPTY,
                        quality=dict(inner=QUALITY, outer=QUALITY,
                                     outer_removed_pct=0))
        else:
            return route.abort()
        route.fulfill(body=json.dumps(data), content_type='application/json')

    def assert_numbers_equal(self, before, after, path=''):
        if isinstance(before, dict):
            self.assertEqual(before.keys(), after.keys(), path)
            for key in before:
                self.assert_numbers_equal(before[key], after[key], f'{path}.{key}')
        elif isinstance(before, list):
            self.assertEqual(len(before), len(after), path)
            for i, value in enumerate(before):
                self.assert_numbers_equal(value, after[i], f'{path}[{i}]')
        elif isinstance(before, (int, float)):
            self.assertAlmostEqual(before, after, delta=1e-6, msg=path)
        else:
            self.assertEqual(before, after, path)

    def assert_view_equal(self, before, after):
        for id in ['inner', 'outer']:
            for key in ['camera', 'ranges', 'aspect']:
                self.assert_numbers_equal(before[id][key], after[id][key],
                                          f'{id}.{key}')

    def test_sliders_preserve_rotated_linked_views(self):
        page = self.browser.new_page(viewport=dict(width=1600, height=1100))
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.route('http://cardiac-cut.test/**', self.route)
        try:
            page.goto('http://cardiac-cut.test/')
            page.wait_for_function("lastPreview && byId('status').textContent.startsWith('Preview ready')")
            for mode in ['orbit', 'turntable']:
                title = 'Orbital rotation' if mode == 'orbit' else 'Turntable rotation'
                page.locator(f'#inner [data-title="{title}"]').click(force=True)
                page.wait_for_function('cameraSyncBusy === 0')
                for source in ['inner', 'outer']:
                    for slider in ['off', 'az', 'el']:
                        with self.subTest(mode=mode, source=source, slider=slider):
                            page.locator('[data-view="fit"]').click()
                            page.wait_for_function('cameraSyncBusy === 0')
                            # Let the camera controller finish its preset transition.
                            page.wait_for_timeout(500)
                            before_drag = page.evaluate(SNAPSHOT)
                            box = page.locator('#' + source).bounding_box()
                            x = box['x'] + box['width'] * .5
                            y = box['y'] + box['height'] * .5
                            page.mouse.move(x, y)
                            page.mouse.down()
                            # Turntable consumes pointer motion on animation frames.
                            for step in range(1, 7):
                                page.mouse.move(x + 70 * step / 6, y + 35 * step / 6)
                                page.wait_for_timeout(40)
                            page.mouse.up()
                            try:
                                page.wait_for_function("""() => {
                                    if (cameraSyncBusy || cameraFollowFrame !== null) return false;
                                    const a = byId('inner')._fullLayout.scene._scene.getCamera();
                                    const b = byId('outer')._fullLayout.scene._scene.getCamera();
                                    return ['eye','up','center'].every(key =>
                                        ['x','y','z'].every(axis => Math.abs(a[key][axis]-b[key][axis]) < 1e-6));
                                }""", timeout=3000)
                            except BrowserTimeout:
                                self.fail(f'Linked cameras did not settle: {page.evaluate(SNAPSHOT)}')
                            before = page.evaluate(SNAPSHOT)
                            self.assertNotEqual(before_drag[source]['camera'],
                                                before[source]['camera'])
                            self.assert_numbers_equal(before['inner']['camera'],
                                                      before['outer']['camera'], 'sync')
                            with page.expect_response(lambda r: '/api/preview' in r.url):
                                page.evaluate("""id => {
                                    const input = byId(id);
                                    input.value = Number(input.value) + Number(input.step);
                                    input.dispatchEvent(new Event('input', {bubbles:true}));
                                }""", slider)
                                immediate = page.evaluate(SNAPSHOT)
                            page.wait_for_function("byId('status').textContent.startsWith('Preview ready')")
                            after = page.evaluate(SNAPSHOT)
                            self.assert_view_equal(before, immediate)
                            self.assert_view_equal(before, after)
                            self.assertNotEqual(before['outer']['plane'], after['outer']['plane'])
                            self.assertNotEqual(before['inner']['mesh'], after['inner']['mesh'])
            self.assertEqual(errors, [])
        finally:
            page.close()


if __name__ == '__main__':
    unittest.main()
