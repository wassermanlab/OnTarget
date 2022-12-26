import React from 'react';
import Select from 'react-select'

function SelectRegion(props) {
    // props
    // page
    // regionType
    // plusMinusGene
    // customCoordinateStart
    // customCoordinateEnd
    // hg19Chrom
    // mm10Chrom
    // hg19Genes
    // mm10Genes
    // geneName

    let chroms = <option></option>
    if (props.genome === "hg19" && props.loadedResources) {
        chroms = Object.keys(props.hg19Chrom).map((item, index) => (
            <option key={index}>{item}</option>
        ))
    } else if (props.genome === "mm10" && props.loadedResources) {
        chroms = Object.keys(props.mm10Chrom).map((item, index) => (
            <option key={index}>{item}</option>
        ))
    }


    let uploadedFiles = <div></div>
    if (props.uploadedFiles !== null) {
        uploadedFiles = <div className='row'>
            <h3>Uploaded Files</h3>
            {props.uploadedFiles.map((item, index) => (
                <div key={index}> {item} </div>))}
        </div>
    }

    return (
        <div>
            <div className='row-margin'>
                <h3>Select Region</h3>
                <p>Note: the maximum size for a region is be 200 kB.</p>
                <div className='row-margin'>
                    <p>Genome</p>
                    <div className="form-check form-check-inline row-margin">
                        <input className="form-check-input" type="radio" name="inlineRadioGenomeOptions" id="inlineRadioGenome1" value="mm10" onChange={props.handleGenomeChange} />
                        <label className="form-check-label" htmlFor="inlineRadioGenome1">mm10</label>
                    </div>
                    <div className="form-check form-check-inline">
                        <input className="form-check-input" type="radio" name="inlineRadioGenomeOptions" id="inlineRadioGenome2" value="hg19" onChange={props.handleGenomeChange} />
                        <label className="form-check-label" htmlFor="inlineRadioGenome2">hg19</label>
                    </div>
                </div>
                {props.genome === "mm10" &&
                    <div>
                        <p>Would you like to liftover to hg19?</p>
                        <div className="form-check form-check-inline row-margin">
                            <input className="form-check-input" type="radio" name="inlineRadioLiftoverOptions" id="inlineRadioLiftoverOptions1" value="true" onChange={props.handleLiftoverChange} />
                            <label className="form-check-label" htmlFor="inlineRadioLiftoverOptions1">Yes</label>
                        </div>
                        <div className="form-check form-check-inline">
                            <input className="form-check-input" type="radio" name="inlineRadioLiftoverOptions" id="inlineRadioLiftoverOptions1" value="false" onChange={props.handleLiftoverChange} />
                            <label className="form-check-label" htmlFor="inlineRadioLiftoverOptions1">No</label>
                        </div>
                    </div>
                }
                <p>Type of region</p>
                <div className="form-check form-check-inline row-margin">
                    <input className="form-check-input" type="radio" name="inlineRadioOptions" id="inlineRadio1" value="geneToGene" onChange={props.handleRegionChange} />
                    <label className="form-check-label" htmlFor="inlineRadio1">Gene to Gene</label>
                </div>
                <div className="form-check form-check-inline">
                    <input className="form-check-input" type="radio" name="inlineRadioOptions" id="inlineRadio2" value="plusMinusBP" onChange={props.handleRegionChange} />
                    <label className="form-check-label" htmlFor="inlineRadio2">+/- n kB from Gene </label>
                </div>
                <div className="form-check form-check-inline">
                    <input className="form-check-input" type="radio" name="inlineRadioOptions" id="inlineRadio3" value="customCoordinates" onChange={props.handleRegionChange} />
                    <label className="form-check-label" htmlFor="inlineRadio3">Custom coordinates around gene</label>
                </div>

                {/* gene name */}
                {(props.regionType === "plusMinusBP" || props.regionType === "geneToGene") && props.genome === "hg19" &&
                    <div className="form-group row-margin">
                        <label htmlFor="selectGene">Gene name</label>
                        <Select onChange={props.handleSelectChange} name="geneName"
                            className="basic-single" classNamePrefix="select" options={props.hg19Genes} />
                    </div>
                }
                {(props.regionType === "plusMinusBP" || props.regionType === "geneToGene") && props.genome === "mm10" &&
                    <div className="form-group row-margin">
                        <label htmlFor="selectGene">Gene name</label>
                        <Select onChange={props.handleSelectChange} name="geneName"
                            className="basic-single" classNamePrefix="select" options={props.mm10Genes} />
                    </div>
                }
                {/* TODO check with oriol what the max is */}
                {props.regionType === "plusMinusBP" &&
                    <div className="form-group row-margin">
                        <label htmlFor="inputBP">n kb from gene</label>
                        <input type="number" className="form-control" min="0" max="100" name="plusMinusGene" placeholder="0"
                            value={props.plusMinusGene} onChange={props.handleChange} />
                    </div>
                }
                {props.regionType === "customCoordinates" &&
                    <div>
                        <div className="form-group row-margin">
                            <label >Chromosome</label>
                            <select onChange={props.handleChange} name="chromosome" value={props.chromosome} className="form-control" id="chromosome">
                                {chroms}
                            </select>
                        </div>
                        <div className="form-group row-margin">
                            <label htmlFor="inputStart">Custom start coordinate (0 based)</label>
                            <input type="number" className="form-control" min="0" name="customCoordinateStart" placeholder="0"
                                value={props.customCoordinateStart} onChange={props.handleChange} />
                        </div>
                        <div className="form-group row-margin">
                            <label htmlFor="inputEnd">Custom end coordinate</label>
                            <input type="number" className="form-control" min="0" name="customCoordinateEnd" placeholder="0"
                                value={props.customCoordinateEnd} onChange={props.handleChange} />
                        </div>
                    </div>
                }
            </div>
            <div className="row-margin">
                <form onSubmit={props.onSubmit}>
                    <h3>Upload Evidence</h3>
                    <p>Only gzipped bed files accepted</p>
                    <div className="form-group">
                        <input type="file" name="evidenceList" onChange={props.onFileChange} multiple />
                    </div>
                    <div className="form-group">
                        <button className="btn btn-primary ontarget-button" type="submit">Upload</button>
                    </div>
                </form>

                {uploadedFiles}

                {/* requestCode */}
                {/* uploadedFiles */}

            </div>
        </div>

    )
}

export default SelectRegion;
