import React from 'react';
import Errors from './errors';
import NavButtons from './navButtons';
import SelectRegion from './selectRegion';
import Loading from './loading';
import axios from 'axios';

class Design extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      errors: [],
      loadedResources: false,
      genes: [],
      enzymes: [],
      tfs: [],
      page: 1,
      //gene
      geneName: "",
      //define_region
      regionType: "",
      genome: "",
      plusMinusGene: 0,
      customCoordinateStart: 0,
      customCoordinateEnd: 0,
      liftover: false,
      evidenceList: '',
      //request_code
      requestCode: "",
      uploadedFiles: null,
      //download, 
      error: null
    };

    {/* Create Mini Promoters */}
    this.onSubmit = this.onSubmit.bind(this);
    this.onFileChange = this.onFileChange.bind(this);
    this.next = this.next.bind(this);
    this.back = this.back.bind(this);
    this.handleRegionChange = this.handleRegionChange.bind(this);
    this.handleChange = this.handleChange.bind(this);
    this.getRegulatoryRegions = this.getRegulatoryRegions.bind(this);
    this.handleGenomeChange = this.handleGenomeChange.bind(this);
    this.handleLiftoverChange = this.handleLiftoverChange.bind(this);
  };

  //set genes, TFs, restrictionEnzymes
  componentDidMount() {
    fetch("http://127.0.0.1:5000/genes_enzymes_tfs") // TODO: change this address
      .then(res => res.json())
      .then(
        (result) => {
          this.setState({
            loadedResources: true,
            genes: result.genes,
            enzymes: result.enzymes,
            tfs: result.tfs
          },
            () => { console.log(this.state); });
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
    console.log("getRegs")
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
      this.setState({requestCode: res.data.request_code, 
      uploadedFiles: res.data.uploaded_files})
    })
  }
  //handles region radio change
  handleChange(event) {
    const target = event.target;
    const value = target.type === 'checkbox' ? target.checked : target.value;
    const name = target.name;

    this.setState({
      [name]: value
    }, () => { console.log(this.state); });
  }

  //handles region radio change
  handleRegionChange(event) {
    this.setState({ regionType: event.target.value }, () => { console.log(this.state); });
  }

  check_errors() {
    // props
    // page
    // regionType
    // plusMinusGene
    // customCoordinateStart
    // customCoordinateEnd
    // genes
    //geneName
    let errors = [];

    if (this.state.page === 1) {
      if (this.state.regionType === "") {
        errors.push("User must select type of region")
      }
      if ((this.state.regionType === "geneToGene" || this.state.regionType === "plusMinusBP") && this.state.geneName === "") {
        errors.push("User must select gene")
      }
      if (this.state.regionType === "plusMinusBP" && (parseInt(this.state.plusMinusGene) > 100 || parseInt(this.state.plusMinusGene) <= 1)) {
        errors.push("n kb from gene must be between 0 and 100")
      }
      //TODO: add error check on custom coordinates
      if (this.state.genome === "") {
        errors.push("User must select genome")
      }
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
    this.setState({ page: currentPage, errors: [] }, () => { console.log(this.state.page); });
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
        <div className="row">
          <div className="col">
          </div>
          <div className="col-8">
              {/* Select Gene */}
              <Loading loadedResources={this.state.loadedResources}></Loading>
              {/* Select Regions */}
              <SelectRegion page={this.state.page}
                handleGenomeChange={this.handleGenomeChange}
                genes={this.state.genes}
                genome={this.state.genome}
                geneName={this.state.geneName}
                handleRegionChange={this.handleRegionChange}
                handleLiftoverChange={this.handleLiftoverChange}
                handleChange={this.handleChange}
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
              {/* buttons */}
              <NavButtons page={this.state.page} next={this.next} back={this.back} getRegulatoryRegions={this.getRegulatoryRegions} />
          </div>
          <div className="col">
          </div>
        </div>
      </div>
    );
  };

}

export default Design