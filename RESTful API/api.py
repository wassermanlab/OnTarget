from flask import Flask
from flask_restful import reqparse, Api, Resource, abort
from organsToSamples import get_sample_names_and_uids_from_organs

app = Flask(__name__)
api = Api(app)

class SamplesFromOrgans(Resource):
    def get(self):
        args = sample_parser.parse_args()
        samples_from_organs = get_sample_names_and_uids_from_organs(args)

        if samples_from_organs is not None:
            return samples_from_organs
        else:
            abort(500)

api.add_resource(SamplesFromOrgans, '/samples')

def _setup_sample_parser():
    global sample_parser
    sample_parser = reqparse.RequestParser()
    sample_parser.add_argument('organs', required=True, type=str, action='append')
    sample_parser.add_argument('X', type=int)
    sample_parser.add_argument('Y', type=int)
    sample_parser.add_argument('cancer', type=bool)
    sample_parser.add_argument('cell_line', type=bool)
    sample_parser.add_argument('treatment', type=bool)

if __name__ == '__main__':
    _setup_sample_parser()
    app.run(debug=True)

