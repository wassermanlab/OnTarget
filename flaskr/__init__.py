import os
import traceback

import requests
from flask import Flask, request, send_from_directory

from flaskr.Endpoints.gene2region import gene2region
from flaskr.Endpoints.gene2sample import gene2sample
from flaskr.Endpoints.name2sample import name2sample
from flaskr.Endpoints.region2enhancer import region2enhancer

remote_host = 'http://gud.cmmt.ubc.ca:8080'


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
    @app.route('/samples/<database>', methods=['GET'])
    def samples(database):
        return name2sample(remote_host, database, request.args.get('name'))

    @app.route('/regions/<database>/<mode>', methods=['POST'])
    def region(database, mode):
        return gene2region(remote_host, database, request.args.get('gene').split(','), request.get_json(), mode)

    @app.route('/gene', methods=['POST'])
    def gene():
        json = request.get_json()
        return gene2sample(samples=json.get("samples", None), cell_line=json.get('cell_line', None),
                           treatment=json.get('treatment', None),
                           XChrom=json.get('XChrom', None), YChrom=json.get('YChrom', None),
                           cancer=json.get('cancer', None))

    @app.route('/enhancer/<database>', methods=['POST'])
    def enhancer(database):
        json = request.get_json()
        # Note I dont know what the default values for some of these should be
        # this currently returns the path of the temp folder created with the data
        return region2enhancer(remote_host, database=database, region=json.get("region", None),
                               length=json.get("length", None), feats=json.get("feats", None),
                               mask_exons=json.get("mask_exons", None), mask_repeats=json.get("mask_repeats", None),
                               ubiquitous=json.get("ubiquitous", None), samples=json.get("samples", []))

    @app.route('/retrieve/<PATH>/<File>', methods=['GET'])
    def retrieve(PATH, File):
        # Current file name options are {chrom}:{start}-{end}.enhancers.bed,
        # {chrom}:{start}-{end}.enhancers.fa, and {chrom}:{start}-{end}.matrix.csv
        send_from_directory(PATH, File, as_attachment=True)

    return app
