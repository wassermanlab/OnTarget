from flask import Flask, flash, request, redirect, url_for, jsonify, send_file
from flask_restful import Resource, Api
from flask_cors import CORS, cross_origin # figure this out on the production server
from werkzeug.utils import secure_filename
from ontarget import hg19, mm10, rest_enzymes, TFs, OnTargetUtils
from ontarget.regions2minips import get_minipromoters, get_minipromoter
from gene2interval import get_intervals_limit_by_gene, get_intervals_limit_by_distance
from interval2regions import get_regions
import os, string, random, distutils
from distutils import util
import shutil
from datetime import datetime




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
def get_minipromoters_request():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        json                = request.json
        requestCode         = str(json["requestCode"]) 
        size                = int(json["size"])                        # put this default on the server side
        selectedEnzymes     = set(json["selectedEnzymes"])
        selectedTFs         = set(json["selectedTFs"])
        regions             = json["regions"]
        designs             = int(json["designs"])                  # put this default on the server side

        # get minipromoters 
        minipromoters = get_minipromoters(regions, designs, selectedEnzymes, size, selectedTFs)

        # make folder for these minimaps 
        currentDateAndTime= datetime.now().strftime("%Y%m%d%H%M%S")
        folder = "minipromoters" + currentDateAndTime
        path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode, folder)
        os.mkdir(path)
        # save the fastas 
        for minip in minipromoters:
            fasta_file = os.path.join(path,
                                      "%s.fa" % "+".join(minip["ids"]))
            OnTargetUtils.write_fasta([minip], fasta_file)
        # download the folder 
        shutil.make_archive(path, 'zip', path)
        zippath = os.path.join(app.config['UPLOAD_FOLDER'], requestCode, (folder +".zip"))

        return send_file(zippath)
    else:
        return {"error":'Content-Type not supported!'}

@app.route('/getminipromoter', methods=["POST"])
def get_minipromoter_request():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        json = request.json
        promoter = json["promoter"]
        enhancers = json["enhancers"]
        requestCode         = str(json["requestCode"])
        minipromoter = get_minipromoter(promoter, enhancers)
        minipromoter = minipromoter.serialize()

        currentDateAndTime= datetime.now().strftime("%Y%m%d%H%M%S")
        filename = "minipromoter" + currentDateAndTime
        path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode)
        # save the fasta
        fasta_file = os.path.join(path,filename+".fa")
        OnTargetUtils.write_fasta([minipromoter], fasta_file)

        return send_file(fasta_file)
    else:
        return {"error":'Content-Type not supported!'}

@app.route('/getregions', methods=["GET"])
def get_regions_request():
    #get args
    regionType          = request.args.get('regionType')            #geneToGene,plusMinusBP,customCoordinates
    genome              = request.args.get('genome')                #hg19 or mm10
    geneName            = request.args.get('geneName')              #name of gene if geneToGene or plusMinusBP
    plusMinusGene       = float(request.args.get('plusMinusGene'))    #how many KB away from gene for plusMinusBP 
    chrom               = request.args.get('chromosome')            #chrom 
    start               = int(request.args.get('customCoordinateStart')) #custom start but then is used for regular start 
    end                 = int(request.args.get('customCoordinateEnd')) #custom end but then is used for regular end 
    requestCode         = request.args.get('requestCode')           #code to directory with evidence
    region_length       = float(request.args.get('region_length'))     #either default or some custom length
    region_score        = float(request.args.get('region_score'))      #either default or some custom score
    cons_score          = float(request.args.get('cons_score'))        #either default or some custom score
    cons_length         = float(request.args.get('cons_length'))       #either default or some custom length
    use_conservation    = bool(distutils.util.strtobool(request.args.get('use_conservation'))) #true or false
    mask_exons          = bool(distutils.util.strtobool(request.args.get('mask_exons')))       #true or false
    mask_repeats        = bool(distutils.util.strtobool(request.args.get('mask_repeats')))     #true or false
    liftover            = request.args.get('liftover') #true or nothing
    
    #fix liftover
    if liftover == "true":
        liftover="hg19"
    #make session
    session =  OnTargetUtils.get_gud_session(genome)
    # get interval if geneToGene or plusMinusBP
    if regionType=="geneToGene":
        interval = get_intervals_limit_by_gene(session, geneName) 
        chrom=interval[0]["chrom"]
        start=interval[0]["start"]
        end=interval[0]["end"]
    elif regionType=="plusMinusBP":
        interval = get_intervals_limit_by_distance(session, geneName, plusMinusGene)
        chrom=interval[0]["chrom"]
        start=interval[0]["start"]
        end=interval[0]["end"]
    # get list of evidence from request code 
    path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode)
    uploadedfiles=os.listdir(path)
    evidence = []
    for i in uploadedfiles:
        evidence.append([os.path.join(app.config['UPLOAD_FOLDER'], requestCode, i), 1.0])
    #get regions
    regions = get_regions(session, chrom, start, end, genome, evidence,
                liftover, region_length,
                region_score,
                use_conservation,
                cons_score,
                cons_length,
                mask_exons, mask_repeats)
    # Write
    bed_file = os.path.join(path, "regions.bed")
    OnTargetUtils.write_bed(regions, bed_file)
    session.close()
    return jsonify(regions)
    
@app.route('/uploadevidence', methods=['POST'])
def upload_file():
    # TODO: check that this is legit! need to accept all bed files under 4MB only save those files that are readable 
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



#this downloads the whole session
@app.route('/download_session', methods=['POST'])
def download_session():
    sessioncode = request.args.get('sessioncode')
    
    path = os.path.join(app.config['UPLOAD_FOLDER'], sessioncode)
    shutil.make_archive(path, 'zip', path)
    zippath = os.path.join(app.config['UPLOAD_FOLDER'], (sessioncode +".zip"))

    return send_file(zippath)


if __name__ == '__main__':
    app.run(debug=True) #remove debug mode for production