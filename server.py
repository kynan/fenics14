"""Launch a livereload server serving up the presentation.

Requires livereload_ (or falls back to SimpleHTTPServer) ::

  pip install git+https://github.com/lepture/python-livereload

.. _livereload: https://github.com/lepture/python-livereload"""

try:
    from livereload import Server
    Server().serve()
except ImportError:
    import SimpleHTTPServer
    import SocketServer

    PORT = 8000
    Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
    httpd = SocketServer.TCPServer(("", PORT), Handler)
    httpd.serve_forever()
