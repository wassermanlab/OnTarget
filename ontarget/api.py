from flask import Flask, flash, request, redirect, url_for
from flask_restful import Resource, Api
from flask_cors import CORS, cross_origin # figure this out on the production server
import os
from werkzeug.utils import secure_filename

UPLOAD_FOLDER = '/home/tamar/Desktop/OnTarget/data/uploads' #change this 

app = Flask(__name__)
CORS(app) # figure out the CORS
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'abcd2435'
# app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
api = Api(app)

def allowed_file(filename):
    print(filename)
    return '.' in filename and filename.endswith(".bed.gz")

@app.route('/uploadevidence', methods=['POST'])
def upload_file():
    print(request.files)
    if request.method == 'POST':
        if 'files' not in request.files:
            flash('No file part')
            return {"message": "failed"}
        files = request.files.getlist('files')
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        return {"message": "passed"}


@app.route('/genes_enzymes_tfs', methods=['GET'])
def get_genes_enzymes_tfs():
    return {"genes": ["gene1", "gene2", "gene3"], 
        "tfs": ["tf1", "tf2", "tf3"],
         "enzymes": ["enzyme1", "enzyme2", "enzyme3"]}

if __name__ == '__main__':
    app.run(debug=True) #remove debug mode for production