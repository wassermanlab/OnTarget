from flask import Flask, flash, request, redirect, url_for, jsonify
from flask_restful import Resource, Api
from flask_cors import CORS, cross_origin # figure this out on the production server
import os
from werkzeug.utils import secure_filename
import string
import random
from ontarget import hg19, mm10, rest_enzymes, TFs
from ontarget.regions2minips import get_minipromoters, get_minipromoter

UPLOAD_FOLDER = '/home/tamar/Desktop/OnTarget/data/uploads' #change this 

app = Flask(__name__)
CORS(app) # figure out the CORS
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'abcd2435' # TODO change this 
# app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
api = Api(app)

def allowed_file(filename):
    print(filename)
    return '.' in filename and filename.endswith(".bed.gz")


def id_generator(size=6, chars=string.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))

@app.route('/getminipromoters', methods=["POST"])
def get_mini_promoters():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        json = request.json
        size = int(json["size"])
        selectedEnzymes = set(json["selectedEnzymes"])
        selectedTFs = set(json["selectedTFs"])
        regions = json["regions"]
        minipromoters = get_minipromoters(regions, size=size, enzymes=selectedEnzymes, tfs= selectedTFs)
        return jsonify(minipromoters)
    else:
        return {"error":'Content-Type not supported!'}

@app.route('/getminipromoter', methods=["POST"])
def get_mini_promoter():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        json = request.json
        promoter = json["promoter"]
        enhancers = json["enhancers"]
        minipromoter = get_minipromoter(promoter, enhancers)
        return jsonify(minipromoter)
    else:
        return {"error":'Content-Type not supported!'}

@app.route('/getregions', methods=["POST"])
def get_regions():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        json = request.json
        ## keys in json for this 
        # genome
        # liftover
        # plusMinusGene
        # customCoordinateStart
        # customCoordinateEnd
        # requestCode
        # regionType
        # geneName

        return {"TODO": "TODO"}
    else:
        return {"error":'Content-Type not supported!'}

@app.route('/uploadevidence', methods=['POST'])
def upload_file():
    # TODO: check that this is legit!
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


@app.route('/genes', methods=['GET'])
def get_genes():
    return {"hg19Chroms": hg19["chroms"],
            "mm10Chroms": mm10["chroms"],
            "hg19Genes": list(hg19["genes"]),
            "mm10Genes": list(mm10["genes"])}


# hg19, mm10, rest_enzymes, TFs
@app.route('/enzymes_tfs', methods=['GET'])
def get_enzymes_tfs():
    tf = list(TFs)
    enzymes = list(rest_enzymes)
    return {"tfs": tf,
         "enzymes": enzymes}

if __name__ == '__main__':
    app.run(debug=True) #remove debug mode for production