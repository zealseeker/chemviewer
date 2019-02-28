from chemviewer.core import app
from chemviewer.core import index

@app.route('/')
def entrence():
    return index('webserver')


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
