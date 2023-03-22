import React from 'react';
import RegionList from './regionList';
import SelectedRegionList from './selectedRegionList';
import Errors from './errors';
import Select from 'react-select'
import axios from 'axios';
import { host } from './host'


class GetRegions extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            designErrors: [],
            loadedResources: false,
            enzymes: [],
            tfs: [],
            //design
            size: 1900, //(in bp; default = 1900; min = 100; max = 4500) 
            designs: 5, // (default = 5; min = 1; max = 100)
            selectedEnzymes: [],
            selectedTFs: [],
            //selected regions to make mini promoter
            selectionError: [],
            selectedRegions: [],
            //download, 
            error: null
        };
        this.selectRegion = this.selectRegion.bind(this);
        this.removeSelectedRegion = this.removeSelectedRegion.bind(this);
        this.clearDesignedPromoter = this.clearDesignedPromoter.bind(this);
        this.downloadDesignedPromoter = this.downloadDesignedPromoter.bind(this);
        this.handleChange = this.handleChange.bind(this);
        this.handleMultiSelectChange = this.handleMultiSelectChange.bind(this);
        this.downloadSuggestedMiniPromoters = this.downloadSuggestedMiniPromoters.bind(this);
    };

    //set genes, TFs, restrictionEnzymes
    componentDidMount() {
        fetch(host+"enzymes_tfs") // TODO: change this address
            .then(res => res.json())
            .then(
                (result) => {
                    let enzymes = []
                    let tfs = []
                    result.enzymes.map((e) => enzymes.push({ value: e, label: e }))
                    result.tfs.map((e) => tfs.push({ value: e, label: e }))
                    this.setState({
                        loadedResources: true,
                        enzymes: enzymes,
                        tfs: tfs
                    });
                },
                (error) => {
                    this.setState({
                        isLoaded: true,
                        error
                    });
                }
            )
    }
    // events
    downloadSuggestedMiniPromoters() {
        let designErrors = []
        if (this.state.size > 4500 || this.state.size < 100) {
            designErrors.push("The maximum size of a MiniPromoter must be between 100-4500bp")
        } if (this.state.designs > 100 || this.state.designs < 1) {
            designErrors.push("The number of requested designs must be between 1-100")
        }
        if (designErrors.length > 0) {
            this.setState({
                designErrors: designErrors
            })
        } else {
            const FileDownload = require('js-file-download');

            const data = JSON.stringify({
                requestCode: this.props.requestCode,
                size: this.state.size,
                designs: this.state.designs,
                selectedEnzymes: this.state.selectedEnzymes,
                selectedTFs: this.state.selectedTFs,
                regions: this.props.regions
            });

            const config = {
                method: 'post',
                url: host+'getminipromoters',  //todo change url
                headers: {
                    'Content-Type': 'application/json'
                },
                data: data,
                responseType: 'blob'
            };

            axios(config)
                .then((response) => {
                    FileDownload(response.data, 'minipromoters.zip');
                });
        }
    };

    selectRegion(index, type) {
        //TODO something
        // get region
        const reg = this.props.regions.filter(r => r.type === type)[index];
        let selectionError = [];
        let selectedRegions = this.state.selectedRegions;
        if (reg.type === "Promoter") {
            if (selectedRegions.filter(r =>
                r.type === "Promoter").length > 0) {
                    selectionError.push("Only one promoter is allowed, remove existing promoter from design before adding a new one.")
                this.setState({ selectionError: selectionError })
                return;
            }
        }
        if (!this.state.selectedRegions.includes(reg)) {
            selectedRegions.push(reg)
            this.setState({
                selectedRegions: selectedRegions,
                selectionError: selectionError
            })
        }

    }
    removeSelectedRegion(index) {
        let selectedRegions = this.state.selectedRegions;
        selectedRegions.splice(index, 1);
        this.setState({ selectedRegions: selectedRegions, selectionError: [] })
    }
    clearDesignedPromoter() {
        this.setState({ selectedRegions: [], selectionError: [] })
    }
    downloadDesignedPromoter() {
        const selectedRegions = this.state.selectedRegions;
        let enhancers = []
        let promoter = 0;
        //sort selected regions
        selectedRegions.forEach((reg) => {
            if (reg.type === "Promoter"){
                promoter = reg;
            } else {
                enhancers.push(reg);
            }
        });
        if (promoter === 0){
            this.setState({
                selectionError: ["MiniPromoter designs must have at least one promoter"]
            });
        } else{
            // TODO: change address
            const FileDownload = require('js-file-download');

            const data = JSON.stringify({
                requestCode: this.props.requestCode,
                enhancers: enhancers,
                promoter: promoter
            });

            const config = {
                method: 'post',
                url: host+'getminipromoter',  //todo change url
                headers: {
                    'Content-Type': 'application/json'
                },
                data: data,
                responseType: 'blob'
            };

            axios(config)
                .then((response) => {
                    FileDownload(response.data, 'minipromoter.fa');
                }).then(()=>{
                    this.setState({ selectedRegions: [], selectionError: [] })
                });
        }
    }
    handleChange(event) {
        const target = event.target;
        const value = target.type === 'checkbox' ? target.checked : target.value;
        const name = target.name;

        this.setState({ [name]: value });
    }
    handleMultiSelectChange(event, name) {
        const values = event.map((item) => { return item.value });
        this.setState({ [name]: values });
    }

    render() {
        return (
            <div className="container-fluid">
                <div className="row">
                    <div className="col">
                    </div>
                    <div className="col-8">
                        <h3>Get suggested regulatory regions</h3>
                        <p>Instructions: To get suggested minipromoters indicate the maximum size of your mini promoter and
                            how many suggested MiniPromoters you would like
                            to design and click "Download Suggested Mini Promoters". </p>
                        {this.props.requestCode==="example" &&
                            <p>
                                For this example we recommend selecting the TFs: <b>NR4A2 & PITX3</b> and the Enzymes: <b>AscI & FseI</b> 
                            </p>
                        }
                        <div className="row" id="numDesigns">
                            <label className="col-sm-4 col-form-label">Number of MiniPromoter Designs (1-100)</label>
                            <div className="col-sm-4">
                                <input type="number" className="form-control" name="designs"
                                    onChange={this.handleChange} value={this.state.designs} />
                            </div>
                        </div>
                        <div className="row" id="numDesigns">
                            <label className="col-sm-4 col-form-label">Maximum size for MiniPromoter (100-4500bp)</label>
                            <div className="col-sm-4">
                                <input type="number" className="form-control" name="size"
                                    value={this.state.size} onChange={this.handleChange} />
                            </div>
                        </div>
                        <div className="row" id="numDesigns">
                            <label className="col-sm-4 col-form-label">TFs to include</label>
                            <div className="col-sm-4">
                                <Select isMulti onChange={e => this.handleMultiSelectChange(e, "selectedTFs")} id="TFSelector"
                                    options={this.state.tfs} className="basic-multi-select" classNamePrefix="select"></Select>
                            </div>
                        </div>
                        <div className="row" id="numDesigns">
                            <label className="col-sm-4 col-form-label">Restriction Enzymes to avoid</label>
                            <div className="col-sm-4">
                                <Select isMulti onChange={e => this.handleMultiSelectChange(e, "selectedEnzymes")} id="EnzymeSelector"
                                    options={this.state.enzymes} className="basic-multi-select" classNamePrefix="select"></Select>
                            </div>
                        </div>
                        <Errors errors={this.state.designErrors} />
                        <button type="button" onClick={this.downloadSuggestedMiniPromoters} className="btn btn-primary ontarget-button" >Download Suggested Mini Promoters</button>
                        <hr /> {/* SPLIT */}
                        <h3>Design Mini promoter from regulatory regions</h3>
                        <p>Instructions: To create a custom minipromoter select at least one promoter and one or more enhancers and click "Get Mini Promoter". </p>
                        <h4>Promoters</h4>
                        <div className='table-wrapper'>
                            <table className="table table-sm table-responsive">
                                <thead>
                                    <tr>
                                        <th scope="col">Chromosome</th>
                                        <th scope="col">Strand</th>
                                        <th scope="col">Start</th>
                                        <th scope="col">End</th>
                                        <th scope="col">Score</th>
                                        <th scope="col">Select Region</th>
                                    </tr>
                                </thead>
                                <RegionList regions={this.props.regions} type="Promoter" selectRegion={this.selectRegion}></RegionList>
                            </table>
                        </div>
                        <h4>Enhancers</h4>
                        <div className='table-wrapper'>
                            <table className="table table-sm table-responsive">
                                <thead>
                                    <tr>
                                        <th scope="col">Chromosome</th>
                                        <th scope="col">Strand</th>
                                        <th scope="col">Start</th>
                                        <th scope="col">End</th>
                                        <th scope="col">Score</th>
                                        <th scope="col">Select Region</th>
                                    </tr>
                                </thead>
                                <RegionList regions={this.props.regions} type="Enhancer" selectRegion={this.selectRegion}></RegionList>
                            </table>
                        </div>
                        <h4>Designed Mini Promoter</h4>
                        <div className='table-wrapper'>
                            <table className="table table-sm table-responsive">
                                <thead>
                                    <tr>
                                        <th scope="col">Region Type</th>
                                        <th scope="col">Chromosome</th>
                                        <th scope="col">Strand</th>
                                        <th scope="col">Start</th>
                                        <th scope="col">End</th>
                                        <th scope="col">Score</th>
                                        <th scope="col">Select Region</th>
                                    </tr>
                                </thead>
                                <SelectedRegionList selectedRegions={this.state.selectedRegions} removeSelectedRegion={this.removeSelectedRegion}></SelectedRegionList>
                            </table>
                        </div>


                        <Errors errors={this.state.selectionError} />
                        <button type="button" onClick={this.downloadDesignedPromoter} className="btn btn-primary ontarget-button" >Download Mini Promoter</button>
                        <button type="button" onClick={this.clearDesignedPromoter} className="btn btn-primary ontarget-button" >Clear Mini Promoter</button>



                    </div>
                    <div className="col">
                    </div>
                </div>
            </div>
        );
    };

}

export default GetRegions