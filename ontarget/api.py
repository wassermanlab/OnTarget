from flask import Flask, flash, request, redirect, url_for
from flask_restful import Resource, Api
from flask_cors import CORS, cross_origin # figure this out on the production server
import os
from werkzeug.utils import secure_filename
import string
import random

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


def id_generator(size=6, chars=string.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))


@app.route('/uploadevidence', methods=['POST'])
def upload_file():
    print(request.files)
    if request.method == 'POST':
        if 'files' not in request.files:
            return {"message": "failed"}
        # get existing directories
        existing_dir = os.listdir(app.config['UPLOAD_FOLDER'])
        # make directory code
        code = ''.join(random.choice(string.ascii_uppercase) for _ in range(6))
        # get random code that doesnt exist
        while (code in existing_dir): 
            code = ''.join(random.choice(string.ascii_uppercase) for _ in range(6))
        path = os.path.join(app.config['UPLOAD_FOLDER'], code)
        os.mkdir(path)
        files = request.files.getlist('files')
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(path, filename))
        uploadedfiles=os.listdir(path)
        return {"message": "passed", 
                "request_code": code, 
                "uploaded_files": uploadedfiles}


@app.route('/genes_enzymes_tfs', methods=['GET'])
def get_genes_enzymes_tfs():
    return {"genes": ["gene1", "gene2", "gene3"], 
        "tfs": ["tf1", "tf2", "tf3"],
         "enzymes": ["enzyme1", "enzyme2", "enzyme3"]}

if __name__ == '__main__':
    app.run(debug=True) #remove debug mode for production