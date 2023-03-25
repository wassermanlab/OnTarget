import React from 'react';
import axios from 'axios';
import { host } from './host'
import GetRegions from "./getRegions";
import adora2JSON from './adora2.json'

class Design extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            regions: adora2JSON,
            requestCode: "adora2example",
            uploadedFilesList: ["dna_accessibility.bed.gz",
            "histone_modification-H3K27ac.bed.gz",
            "histone_modification-H3K36me3.bed.gz",
            "histone_modification-H3K4me1.bed.gz",
            "histone_modification-H3K4me3.bed.gz",
            "tf_binding-POL2.bed.gz"],
        };

        this.downloadEvidence = this.downloadEvidence.bind(this);
    }

    downloadEvidence() {
        const FileDownload = require('js-file-download');

        const data = JSON.stringify({
            requestCode: this.state.requestCode,
        });

        const config = {
            method: 'post',
            url: host + 'getadora2evidence',
            headers: {
                'Content-Type': 'application/json'
            },
            data: data,
            responseType: 'blob'
        };

        axios(config)
            .then((response) => {
                FileDownload(response.data, 'adora2evidence.zip');
            });
    }

    render() {
        return (
            <div className="container-fluid">

                <div className="row">
                    <div className="col">
                    </div>
                    <div className="col-8">
                        <h3>Select Region</h3>
                        <p>Note: the maximum size for a region is be 200 kB.</p>
                        <div className='row-margin'>
                            <p>Genome</p>
                            <div className="form-check form-check-inline row-margin">
                                <input className="form-check-input" type="radio" name="inlineRadioGenomeOptions" id="inlineRadioGenome1" value="mm10" checked disabled />
                                <label className="form-check-label" htmlFor="inlineRadioGenome1">mm10</label>
                            </div>
                            <div className="form-check form-check-inline">
                                <input className="form-check-input" type="radio" name="inlineRadioGenomeOptions" id="inlineRadioGenome2" value="hg19" disabled />
                                <label className="form-check-label" htmlFor="inlineRadioGenome2">hg19</label>
                            </div>
                        </div>
                        <div>
                            <p>Would you like to liftover to hg19?</p>
                            <div className="form-check form-check-inline row-margin">
                                <input className="form-check-input" type="radio" name="inlineRadioLiftoverOptions" id="inlineRadioLiftoverOptions1" value="true" checked disabled />
                                <label className="form-check-label" htmlFor="inlineRadioLiftoverOptions1">Yes</label>
                            </div>
                            <div className="form-check form-check-inline">
                                <input className="form-check-input" type="radio" name="inlineRadioLiftoverOptions" id="inlineRadioLiftoverOptions1" value="false" disabled />
                                <label className="form-check-label" htmlFor="inlineRadioLiftoverOptions1">No</label>
                            </div>
                        </div>
                        <p>Type of region</p>
                        <div className="form-check form-check-inline row-margin">
                            <input className="form-check-input" type="radio" name="inlineRadioOptions" id="inlineRadio1" value="geneToGene" disabled />
                            <label className="form-check-label" htmlFor="inlineRadio1">Gene to Gene</label>
                        </div>
                        <div className="form-check form-check-inline">
                            <input className="form-check-input" type="radio" name="inlineRadioOptions" id="inlineRadio2" value="plusMinusBP" checked disabled />
                            <label className="form-check-label" htmlFor="inlineRadio2">+/- n kB from Gene </label>
                        </div>
                        <div className="form-check form-check-inline">
                            <input className="form-check-input" type="radio" name="inlineRadioOptions" id="inlineRadio3" value="customCoordinates" disabled />
                            <label className="form-check-label" htmlFor="inlineRadio3">Custom Coordinates</label>
                        </div>
                        <div className="form-group row-margin">
                            <label htmlFor="selectGene">Gene name</label>
                            <input type="text" className="form-control" placeholder="ADORA2" disabled />
                        </div>
                        <div className="form-group row-margin">
                            <label htmlFor="inputBP">n kb from gene</label>
                            <input type="number" className="form-control" defaultValue="100" disabled />
                        </div>
                        <hr />
                        <div className="row-margin">
                            <h3>Upload Evidence</h3>
                            <p>Only .bed.gz or .bed files less than 4 MB accepted. Successfully uploaded files will be displyed after pressing the "upload" button.</p>
                            <div className='row'>
                                <h3>Uploaded Files</h3>
                                {this.state.uploadedFilesList.map((item, index) => (
                                    <div key={index}> {item} </div>))}
                            </div>
                        </div>
                        <button className="btn btn-primary ontarget-button" onClick={this.downloadEvidence}>Download Example Evidence</button>
                        <hr />
                        <h3>Advanced Parameters</h3>
                        <div className="row-margin">
                            <div className="row advanceParam">
                                <label className="col-sm-4 col-form-label">Region length (0-1000)</label>
                                <div className="col-sm-4">
                                    <input type="number" className="form-control" value="100" name="region_length"
                                         disabled />
                                </div>
                            </div>
                            <div className="row advanceParam">
                                <label className="col-sm-4 col-form-label">Region score (0-1)</label>
                                <div className="col-sm-4">
                                    <input type="number" className="form-control" min="0" name="region_score"
                                        value="0.5546155574351114" disabled/>
                                </div>
                            </div>
                            <div className="row advanceParam">
                                <label className="col-sm-4 col-form-label"> Conserved region length (0-1000)</label>
                                <div className="col-sm-4">
                                    <input type="number" className="form-control" min="0" name="cons_length"
                                        value="10" disabled />
                                </div>
                            </div>
                            <div className="row advanceParam">
                                <label className="col-sm-4 col-form-label">Conserved region score (0-1)</label>
                                <div className="col-sm-4">
                                    <input type="number" className="form-control" min="0" name="cons_score"
                                        value="0.6" disabled />
                                </div>
                            </div>
                            <div className="row advanceParam">
                                <label className="col-sm-3 col-form-label">Use conservation</label>
                                <div className="col-sm-8">
                                    <input type="checkbox" className="form-check-input" min="0" name="use_conservation"
                                        checked disabled />
                                </div>
                            </div>
                            <div className="row advanceParam">
                                <label className="col-sm-3 col-form-label">Mask Exons</label>
                                <div className="col-sm-8">
                                    <input type="checkbox" className="form-check-input" min="0" name="mask_exons"
                                        checked disabled />
                                </div>
                            </div>
                            <div className="row advanceParam">
                                <label className="col-sm-3 col-form-label">Mask Repeats</label>
                                <div className="col-sm-8">
                                    <input type="checkbox"  className="form-check-input" min="0" name="mask_repeats"
                                        checked disabled />
                                </div>
                            </div>
                            <hr />
                        </div>
                    </div>
                    <div className="col">
                    </div>
                </div>

                <GetRegions requestCode={this.state.requestCode} regions={this.state.regions}></GetRegions>
            </div>
        );
    };

}

export default Design
