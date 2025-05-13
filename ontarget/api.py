from flask import Flask, flash, request, redirect, url_for, jsonify, send_file
from flask_restful import Resource, Api
from flask_cors import CORS, cross_origin # figure this out on the production server
from werkzeug.utils import secure_filename
from ontarget import hg19, mm10, rest_enzymes, TFs, OnTargetUtils
from ontarget.regions2minips import get_minipromoters, get_minipromoter
from ontarget.gene2interval import get_intervals_limit_by_gene, get_intervals_limit_by_distance
from ontarget.interval2regions import get_regions
import os, string, random, distutils
from distutils import util
import shutil
from datetime import datetime
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from pybedtools import BedTool
import gzip
import re

app = Flask(__name__)
app.config.from_pyfile('config.py')
#limiter = Limiter(
#    app,
#    key_func=get_remote_address,
#    default_limits=["500 per day", "5 per second"]
#)

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

cors = CORS(app, resources={r"/*": {"origins": ["http://localhost:3000",
						"http://gud.cmmt.ubc.ca",
						"http://gud.cmmt.ubc.ca:8080",
						"http://ontarget.cmmt.ubc.ca",
						"http://ontarget.cmmt.ubc.ca:8080"]}})
#app.config["CORS_HEADERS"] = "Content-Type"
api = Api(app)

def allowed_file(filename):
    return '.' in filename and (filename.endswith(".bed.gz") or filename.endswith(".bed"))


def id_generator(size=6, chars=string.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))

@app.route('/getminipromoters', methods=["POST"])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_minipromoters_request():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        # designs : int , max. number of MiniPromoters to design per promoter region (default = 5; min = 1; max = 100)
        # size : int , max. MiniPromoter size (in bp; default = 1900; min = 100; max = 4500) 
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

        #check if its for the example
        path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode)
        pathExists = os.path.exists(path)
        if requestCode == "example" and not pathExists: 
            os.makedirs(path)
        path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode, folder)
        os.makedirs(path)
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
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_minipromoter_request():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        json = request.json
        promoter = json["promoter"]
        enhancers = json["enhancers"]
        requestCode         = str(json["requestCode"])
        minipromoter = get_minipromoter(promoter, enhancers)
        minipromoter = minipromoter.serialize()

        #check if its for the example
        path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode)
        pathExists = os.path.exists(path)
        if requestCode == "example" and not pathExists: 
            os.makedirs(path)
        
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
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_regions_request():
    #get args
    regionType          = request.args.get('regionType')            #geneToGene,plusMinusBP,customCoordinates
    genome              = request.args.get('genome')                #hg19 or mm10
    geneName            = request.args.get('geneName')              #name of gene if geneToGene or plusMinusBP
    plusMinusGene       = int(request.args.get('plusMinusGene'))    #how many KB away from gene for plusMinusBP 
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
    else:
        liftover=None
    #make session
    session =  OnTargetUtils.get_gud_session(genome)
    # get interval if geneToGene or plusMinusBP
    if regionType=="geneToGene":
        interval = get_intervals_limit_by_gene(session, geneName, genome) 
        chrom=interval[0]["chrom"]
        start=interval[0]["start"]
        end=interval[0]["end"]
    elif regionType=="plusMinusBP":
        interval = get_intervals_limit_by_distance(session, geneName, (plusMinusGene * 1000))
        chrom=interval[0]["chrom"]
        start=interval[0]["start"]
        end=interval[0]["end"]
    # get list of evidence from request code 
    if requestCode == "":
        existing_dir = os.listdir(app.config['UPLOAD_FOLDER'])
        # make directory code
        code = ''.join(random.choice(string.ascii_uppercase) for _ in range(6))
        # get random code that doesnt exist
        while (code in existing_dir): 
            code = ''.join(random.choice(string.ascii_uppercase) for _ in range(6))
        path = os.path.join(app.config['UPLOAD_FOLDER'], code)
        os.mkdir(path)
    else:
        path = os.path.join(app.config['UPLOAD_FOLDER'], requestCode)
    uploadedfiles=os.listdir(path)
    evidence = []
    for i in uploadedfiles:
        evidence.append([os.path.join(app.config['UPLOAD_FOLDER'], requestCode, i), 1.0])

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

def get_size(fobj):
    if fobj.content_length:
        return fobj.content_length

    try:
        pos = fobj.tell()
        fobj.seek(0, 2)  #seek to end
        size = fobj.tell()
        fobj.seek(pos)  # back to original position
        return size
    except (AttributeError, IOError):
        pass

    # in-memory file object that doesn't support seeking or tell
    return 0  #assume small enough

@app.route('/getpitx3evidence', methods=["POST"])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_pitx3_evidence():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        zippath = os.path.join(app.config['PITX3_EVIDENCE'], "PITX3.zip") 
        return send_file(zippath)
    else:
        return {"error":'Content-Type not supported!'}
    
@app.route('/getadora2evidence', methods=["POST"])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_adora2_evidence():
    content_type = request.headers.get('Content-Type')
    if (content_type == 'application/json'):
        zippath = os.path.join(app.config['PITX3_EVIDENCE'], "ADORA2.zip") 
        return send_file(zippath)
    else:
        return {"error":'Content-Type not supported!'}

@app.route('/uploadevidence', methods=['POST'])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def upload_file():
    # TODO: check that this is legit! need to accept all bed files under 4MB only save those files that are readable 
    if request.method == 'POST':
        if 'files' not in request.files:
            return {"message": "failed"}
        # get existing directories
        # make directory code if not exists
        if not os.path.exists(app.config['UPLOAD_FOLDER']):
            os.mkdir(app.config['UPLOAD_FOLDER'])
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
            
            if not allowed_file(file.filename):
                raise Exception("not allowed extension: "+file.filename)
            elif get_size(file) > 4 * (1024 * 1024):
                print("bed file too big: "+file.filename)
                raise Exception("file too big: "+file.filename)
            else:
                try:
                    file_contents = file.read().decode("utf-8")
                    file.seek(0)
                except UnicodeDecodeError:
                    if file.filename.endswith(".bed.gz"):
                        try:
                            with gzip.open(file.stream, 'rt', encoding='utf-8') as f:
                                file_contents = f.read()
                                file.seek(0)
                        except Exception as e:
                            raise Exception("Error reading gzipped file: {e}")
                    else:
                        raise Exception("Error reading file: {e}")
                except Exception as e:
                    raise Exception("Error reading file: {e}")
                suspicious_patterns = [
                    r"#!",
                    r"\bscript\b", 
                    r"\beval\b",
                    r"\bexec\b",
                    r"rm\s+-rf",
                    r"\bbash\b",  
                    r"\bimport\b",
                    r"\bpython\b",
                    r"\bshell\b",
                    r"\bssh\b",
                    r"\bscp\b",
                    r"\bwget\b",
                    r"\bcurl\b",
                    r"\bnc\b",
                    r"\btelnet\b",
                    r"\bncat\b",
                    r"\bngrok\b",
                    r"\bexecvp\b",
                    r"\bfork\b",
                ]

                for pattern in suspicious_patterns:
                    if re.search(pattern, file_contents):
                        raise Exception(f"Suspicious content found: '{pattern}'") #abort(400)
                filename = secure_filename(file.filename)
                file.seek(0)
                file.save(os.path.join(path, filename))
                # try to read in file if its not a bed file delete the files
                try:
                    b = BedTool(os.path.join(path, filename))
                    #print("file will be validated: ", filename)
                    #print("count----", b.count())
                    #print("head----", b.head())
                    # if b[0].chrom and b[0].chromStart and b[0].chromEnd:
                    #     print("file is readable as bed: ", filename)
                    # else:
                    #     raise Exception("file not readable as bed")
                except:
                    os.remove(os.path.join(path, filename))
                    print("file not readable as bed: ", filename)
                    raise Exception("file not readable as bed")
        uploadedfiles=os.listdir(path)
        return {"message": "passed", 
                "request_code": code, 
                "uploaded_files": uploadedfiles}


@app.route('/genes', methods=['GET'])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_genes():
    return {"hg19Chroms": hg19["chroms"],
            "mm10Chroms": mm10["chroms"],
            "hg19Genes": list(hg19["genes"]),
            "mm10Genes": list(mm10["genes"])}
    #response = jsonify({"hijo": "perra",
    #                    "hg19Chroms": hg19["chroms"],
    #                    "mm10Chroms": mm10["chroms"],
    #                    "hg19Genes": list(hg19["genes"]),
    #                    mm10Genes": list(mm10["genes"])})
    #response.headers.add("Access-Control-Allow-Origin", "*")
    #return response


# hg19, mm10, rest_enzymes, TFs
@app.route('/enzymes_tfs', methods=['GET'])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def get_enzymes_tfs():
    tf = list(TFs)
    enzymes = list(rest_enzymes)
    return {"tfs": tf,
         "enzymes": enzymes}



#this downloads the whole session
@app.route('/download_session', methods=['POST'])
#@cross_origin(origin="*",headers=["Content-Type", "Authorization"])
def download_session():
    sessioncode = request.args.get('sessioncode')
    
    path = os.path.join(app.config['UPLOAD_FOLDER'], sessioncode)
    shutil.make_archive(path, 'zip', path)
    zippath = os.path.join(app.config['UPLOAD_FOLDER'], (sessioncode +".zip"))

    return send_file(zippath)


if __name__ == '__main__':
    app.run(debug=False) #remove debug mode for production
