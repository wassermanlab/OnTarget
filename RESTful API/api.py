from flask import Flask
from flask_restful import reqparse, Api, Resource, abort
from nameToSample import get_sample_for_name

app = Flask(__name__)
api = Api(app)

parser = reqparse.RequestParser()

class Samples(Resource):
    def get(self, genome, name):
        samples_for_name = get_sample_for_name(genome, name)

        if samples_for_name is not None:
            return get_sample_for_name(genome, name)
        else:
            abort(500)
            

api.add_resource(Samples, '/<genome>/samples/<name>')

if __name__ == '__main__':
    app.run(debug=True)