from chemviewer.core import app
from chemviewer.core import index
import os

@app.route('/')
def entrence():
    zealseeker = False
    if 'CHEMVIEWER_SERVER' in os.environ:
        zealseeker = os.environ['CHEMVIEWER_SERVER'] == 'zealseeker'
    return index('webserver', zealseeker=zealseeker)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
