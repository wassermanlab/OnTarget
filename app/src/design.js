import React from 'react';
import Errors from './errors';
import NavButtons from './navButtons';
import CreateMiniPromoter from './createMiniPromoter';
import AddEvidence from './addEvidence';
import SelectRegion from './selectRegion';
import Loading from './loading';

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
      //upload_files
      fileList: [{
        file: "",
        evidence: ""
      }],
      //create_minipromoter
      regulatoryRegions: [
        {
          "chrom": "chr10",
          "start": 103985952,
          "end": 103986304,
          "type": "TSS",
          "id": "RegulatoryRegion1",
          "score": 0.765,
          "strand": -1,
          "qualifiers": null,
          "selected": false
        },
        {
          "chrom": "chr10",
          "start": 103989659,
          "end": 103989823,
          "type": "Enhancer",
          "id": "RegulatoryRegion2",
          "score": 0.765,
          "strand": -1,
          "qualifiers": null,
          "selected": false
        },
        {
          "chrom": "chr10",
          "start": 103989913,
          "end": 103990269,
          "type": "Enhancer",
          "id": "RegulatoryRegion3",
          "score": 0.765,
          "strand": -1,
          "qualifiers": null,
          "selected": false
        }
      ],
      miniPromoterCart: [],
      //download, 
      error: null  
    };

    this.next = this.next.bind(this)
    this.back = this.back.bind(this)
    this.handleRegionChange = this.handleRegionChange.bind(this)
    this.handleChange = this.handleChange.bind(this);
    this.addEvidenceRow = this.addEvidenceRow.bind(this);
    this.removeEvidenceRow = this.removeEvidenceRow.bind(this);
    this.changeEvidence = this.changeEvidence.bind(this);
    this.handleRegCheckBox = this.handleRegCheckBox.bind(this);
    this.addMiniPromoter = this.addMiniPromoter.bind(this);
    this.removeMiniPromoter = this.removeMiniPromoter.bind(this);
    this.getRegulatoryRegions = this.getRegulatoryRegions.bind(this);
    this.handleGenomeChange = this.handleGenomeChange.bind(this);
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
  handleGenomeChange(event){
    this.setState({ genome: event.target.value });
  }

  //sendRegulatoryRegions
  getRegulatoryRegions(){
   console.log("getRegs")
  }

  //removes minipromoter
  removeMiniPromoter(i){
    let miniPromoterCart = this.state.miniPromoterCart;
    miniPromoterCart.splice(i, 1);
    this.setState({ miniPromoterCart });
  }

  //adds minipromoter
  addMiniPromoter() {
    let regulatoryRegions = this.state.regulatoryRegions;
    let miniPromoterCart = this.state.miniPromoterCart;
    //regulatoryRegionsList, totalScore
    let newMiniPromoter = { regulatoryRegionsList: [], totalScore: 0 };
    regulatoryRegions.forEach(function (region, index) {
      if (region.selected) {
        newMiniPromoter.regulatoryRegionsList.push(region.id);
        newMiniPromoter.totalScore += region.score;
      }
      regulatoryRegions[index]["selected"] = false;
    });
    if (newMiniPromoter.regulatoryRegionsList.length > 0) {
      miniPromoterCart.push(newMiniPromoter);
    }
    this.setState({ regulatoryRegions, miniPromoterCart },
      () => { console.log(this.state); });
  }

  //handles checkbox change in handleRegCheckBox
  handleRegCheckBox(i, e) {
    let regulatoryRegions = this.state.regulatoryRegions;
    regulatoryRegions[i][e.target.name] = e.target.checked;
    this.setState({ regulatoryRegions }, () => { console.log(this.state); });
  }

  //handles evidence change 
  changeEvidence(i, e) {
    let fileList = this.state.fileList;
    fileList[i][e.target.name] = e.target.value;
    this.setState({ fileList }, () => { console.log(this.state); });
  }

  //removes another row of evidence
  removeEvidenceRow(i) {
    let fileList = this.state.fileList;
    fileList.splice(i, 1);
    this.setState({ fileList });
  }

  //adds another row of evidence
  addEvidenceRow() {
    this.setState({
      fileList: [...this.state.fileList, {
        file: "",
        evidence: ""
      }]
    }, () => { console.log(this.state); })
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

    if (this.state.page === 1){
      if(this.state.regionType === "" ) {
        errors.push("User must select type of region")
      }
      if((this.state.regionType === "geneToGene" || this.state.regionType === "plusMinusBP") && this.state.geneName==="") {
        errors.push("User must select gene")
      }
      if(this.state.regionType === "plusMinusBP" && (parseInt(this.state.plusMinusGene)>100 || parseInt( this.state.plusMinusGene)<=1)) {
        errors.push("n kb from gene must be between 0 and 100")
      }
      //TODO: add error check on custom coordinates
      if(this.state.genome === "" ) {
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
            <form>
              {/* Select Gene */}
              <Loading loadedResources={this.state.loadedResources}></Loading>
              {/* Select Regions */}
              <SelectRegion page={this.state.page}
                handleGenomeChange={this.handleGenomeChange}
                genes={this.state.genes}
                geneName={this.state.geneName}
                handleRegionChange={this.handleRegionChange}
                handleChange={this.handleChange}
                regionType={this.state.regionType}
                plusMinusGene={this.state.plusMinusGene}
                customCoordinateStart={this.state.customCoordinateStart}
                customCoordinateEnd={this.state.customCoordinateEnd}
              ></SelectRegion>
              {/* Add Evidence */}
              <AddEvidence page={this.state.page} fileList={this.state.fileList}
                addEvidenceRow={this.addEvidenceRow} removeEvidenceRow={this.removeEvidenceRow} changeEvidence={this.changeEvidence}></AddEvidence>
              {/* Create Mini Promoters */}
              <CreateMiniPromoter page={this.state.page} regulatoryRegions={this.state.regulatoryRegions} miniPromoterCart={this.state.miniPromoterCart}
                handleRegCheckBox={this.handleRegCheckBox} addMiniPromoter={this.addMiniPromoter} removeMiniPromoter={this.removeMiniPromoter}></CreateMiniPromoter>
              {/* error handeling */}
              <Errors errors={this.state.errors} />
              {/* buttons */}
              <NavButtons page={this.state.page} next={this.next} back={this.back} getRegulatoryRegions={this.getRegulatoryRegions}/>
            </form>
          </div>
          <div className="col">
          </div>
        </div>
      </div>
    );
  };

}

export default Design