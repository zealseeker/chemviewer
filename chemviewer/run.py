from chemviewer.standalone import app
import os
from sys import argv
import webbrowser

def main():
    if len(argv) > 1:
        webbrowser.open_new('http://127.0.0.1:5000/?file='+os.path.abspath(argv[1]))
    else:
        webbrowser.open_new('http://127.0.0.1:5000/')
    app.run()
