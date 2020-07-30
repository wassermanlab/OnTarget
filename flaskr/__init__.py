import os

from flask import Flask, request

from flaskr.Endpoints.name2sample import name2sample

remote_host='http://gud.cmmt.ubc.ca:8080'

def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # eg host/samples/hg38/?name=endothelial%20cells
    @app.route('/samples/<database>/')
    def samples(database):
        return name2sample(remote_host, database, request.args.get('name'))

    return app