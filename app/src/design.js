import React from 'react';
import Errors from './errors';
import SelectRegion from './design_components/selectRegion';
import Loading from './design_components/loading';
import axios from 'axios';
import GetRegions from "./getRegions";

class Design extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      errors: [],
      page: 1,
      loadedResources: false,
      hg19Chrom: null,
      mm10Chrom: null,
      hg19Genes: null,
      mm10Genes: null,
      //gene
      geneName: "",
      //define_region
      regionType: "",
      genome: "",
      plusMinusGene: 0,
      chromosome: 1,
      customCoordinateStart: 0,
      customCoordinateEnd: 0,
      liftover: false,
      evidenceList: '',
      //request_code
      requestCode: "",
      uploadedFiles: null,
      //
      //componenet did mount error , 
      error: null
    };

    this.onSubmit = this.onSubmit.bind(this);
    this.onFileChange = this.onFileChange.bind(this);
    this.next = this.next.bind(this);
    this.back = this.back.bind(this);
    this.handleRegionChange = this.handleRegionChange.bind(this);
    this.handleChange = this.handleChange.bind(this);
    this.handleSelectChange = this.handleSelectChange.bind(this);
    this.getRegulatoryRegions = this.getRegulatoryRegions.bind(this);
    this.handleGenomeChange = this.handleGenomeChange.bind(this);
    this.handleLiftoverChange = this.handleLiftoverChange.bind(this);
  };

  //set genes, TFs, restrictionEnzymes
  componentDidMount() {
    fetch("http://127.0.0.1:5000/genes") // TODO: change this address
      .then(res => res.json())
      .then(
        (result) => {
          let hg19Genes = []
          let mm10Genes = []
          result.hg19Genes.map((e) => hg19Genes.push({ value: e, label: e }))
          result.mm10Genes.map((e) => mm10Genes.push({ value: e, label: e }))
          this.setState({
            loadedResources: true,
            hg19Genes: hg19Genes,
            mm10Genes: mm10Genes,
            hg19Chrom: result.hg19Chroms,
            mm10Chrom: result.mm10Chroms,
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

  //handleGenomeChange
  handleGenomeChange(event) {
    this.setState({ genome: event.target.value });
  }

  //handleLiftoverChange
  handleLiftoverChange(event) {
    this.setState({ liftover: event.target.value });
  }
  //sendRegulatoryRegions
  getRegulatoryRegions() {
    const errors = this.check_errors()
    if (errors.length !== 0) {
      this.setState({ errors: errors });
      return null;
    } else {
      this.setState({ errors: errors,
        page:2 
      });
      return null;
    }
    // Args to send 
    // ##   regionType // geneToGene,plusMinusBP,customCoordinates
    // ##   genome
    // ##   plusMinusGene
    // ##   chromosome
    // ##   customCoordinateStart
    // ##   customCoordinateEnd
    // ##   liftover // true,false 
    // ##   requestCode
    // ##   geneName
  }
  onFileChange(e) {
    this.setState({ evidenceList: e.target.files })
  }
  onSubmit(e) {
    e.preventDefault()
    var formData = new FormData();
    for (const key of Object.keys(this.state.evidenceList)) {
      formData.append('files', this.state.evidenceList[key])
    }
    axios.post("http://127.0.0.1:5000/uploadevidence", formData, { // TODO change this 
    }).then(res => {
      this.setState({
        requestCode: res.data.request_code,
        uploadedFiles: res.data.uploaded_files
      })
    })
  }
  // handles gene select change 
  handleSelectChange(event) {
    this.setState({ geneName: event.value })
  }
  //handles region radio change
  handleChange(event) {
    const target = event.target;
    const value = target.type === 'checkbox' ? target.checked : target.value;
    const name = target.name;

    this.setState({
      [name]: value
    });
  }

  //handles region radio change
  handleRegionChange(event) {
    this.setState({ regionType: event.target.value });
  }

  check_errors() {
    let errors = [];

    if (this.state.regionType === "") {
      errors.push("User must select type of region")
    }
    if ((this.state.regionType === "geneToGene" || this.state.regionType === "plusMinusBP") && this.state.geneName === "") {
      errors.push("User must select gene")
    }
    if (this.state.regionType === "plusMinusBP" && (parseInt(this.state.plusMinusGene) > 100 || parseInt(this.state.plusMinusGene) <= 1)) {
      errors.push("n kb from gene must be between 0 and 100")
    }
    if (this.state.regionType === "customCoordinates" && this.state.genome === "") {
      errors.push("User must select genome to verify coordinates")
    } else if (this.state.regionType === "customCoordinates") {
      const chromSize = this.state.genome === "hg19" ?
        this.state.hg19Chrom[this.state.chromosome] :
        this.state.mm10Chrom[this.state.chromosome];
      if (parseInt(this.state.customCoordinateEnd) > parseInt(chromSize)) {
        errors.push("End coordinate must be valid for selected chromosome")
      } 
      if ((parseInt(this.state.customCoordinateEnd) - parseInt(this.state.customCoordinateStart)) > 200000) {
        errors.push("End coordinate cannot be more than 200000bp away from start coordinate")
      }
      if (parseInt(this.state.customCoordinateEnd) <= parseInt(this.state.customCoordinateStart)) {
        errors.push("End coordinate must be larger than start coordinate")
      }
      if (parseInt(this.state.customCoordinateEnd) <= 0 || parseInt(this.state.customCoordinateStart) < 0) {
        errors.push("Coordinates must be positive integers")
      }
    }
    if (this.state.genome === "") {
      errors.push("User must select genome")
    }
    if (this.state.requestCode === "" || (this.state.uploadedFiles.length === 0)) {
      errors.push("User must upload at least 1 bed file of evidence.")
    }
    return errors;
  }
  // sets page to page+1
  next() {
    // todo: add error checking
    const errors = this.check_errors()
    if (errors.length !== 0) {
      this.setState({ errors: errors });
      return null;
    }
    let currentPage = this.state.page;
    currentPage = currentPage + 1;
    this.setState({ page: currentPage, errors: [] });
  }
  // sets page to page-1
  back() {
    let currentPage = this.state.page;
    currentPage = currentPage - 1;
    this.setState({ page: currentPage });
  }
  render() {
    return (
      <div className="container-fluid">
        {this.state.page === 1 &&
          <div className="row">
            <div className="col">
            </div>
            <div className="col-8">
              {/* Select Gene */}
              <Loading loadedResources={this.state.loadedResources}></Loading>
              {/* Select Regions */}
              <SelectRegion page={this.state.page}
                loadedResources={this.state.loadedResources}
                handleGenomeChange={this.handleGenomeChange}
                hg19Chrom={this.state.hg19Chrom}
                mm10Chrom={this.state.mm10Chrom}
                hg19Genes={this.state.hg19Genes}
                mm10Genes={this.state.mm10Genes}
                chromosome={this.state.chromsome}
                genome={this.state.genome}
                geneName={this.state.geneName}
                handleRegionChange={this.handleRegionChange}
                handleLiftoverChange={this.handleLiftoverChange}
                handleChange={this.handleChange}
                handleSelectChange={this.handleSelectChange}
                regionType={this.state.regionType}
                liftover={this.state.liftover}
                plusMinusGene={this.state.plusMinusGene}
                customCoordinateStart={this.state.customCoordinateStart}
                customCoordinateEnd={this.state.customCoordinateEnd}
                onSubmit={this.onSubmit} onFileChange={this.onFileChange}
                evidenceList={this.state.evidenceList}
                requestCode={this.state.requestCode}
                uploadedFiles={this.state.uploadedFiles}
              ></SelectRegion>
              {/* error handeling */}
              <Errors errors={this.state.errors} />
              <button to="/getregions" onClick={this.getRegulatoryRegions} className="btn btn-primary ontarget-button">Get Regulatory Regions</button>

            </div>
            <div className="col">
            </div>
          </div>}
        {this.state.page === 2 &&
          <GetRegions requestCode={this.state.requestCode}></GetRegions>}
      </div>
    );
  };

}

export default Design