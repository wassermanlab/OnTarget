import React, { useState, useEffect } from 'react';
// import IGV from './igv'
import { ReactSortable } from "react-sortablejs";
import './design.css';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome'
import { faSyncAlt } from '@fortawesome/free-solid-svg-icons'

function Construct(props) {
  // regulatoryRegions, genomeSize
  const [c, setC] = useState([])


  useEffect(() => {
    let offset = 0;
    let newC = []
    props.regulatoryRegions.forEach(e => {
      if (e.inConstruct) {
        const lenPercentage = (((e.end + e.after_end) - (e.start - e.before_start)) / props.genomeSize) * 100;
        const y2 = offset + lenPercentage;
        newC.push({ "type": e.type, "name": e.name, "y1": offset, "y2": y2 })
        offset = y2;
      }
    })
    setC(newC);
  }, [props.regulatoryRegions, props.genomeSize]);
  let offset = 0;
  return (
    <g>
      {c.map((item, index) => {
        return (
          <line key={index} x1={`${item.y1}%`} y1="20" x2={`${item.y2}%`}
            y2="20" style={{ stroke: (item.type === "promoter" ? "#2d7863" : "#6bcaaf"), strokeWidth: 10 }} />)
      })}
    </g>)
}


export default function Design(props) {
  const [genomeSize, setGenomeSize] = useState(4500);
  const [cargoSize, setCargoSize] = useState(2500);
  const [designSize, setDesignSize] = useState(2000);
  const [regulatoryRegions, setRegulatoryRegions] = useState([
    { "name": "p1", "type": "promoter", "score": .7, "start": 1, "end": 100, "before_start": 0, "after_end": 0, "inConstruct": false },
    { "name": "p2", "type": "promoter", "score": .5, "start": 200, "end": 300, "before_start": 0, "after_end": 0, "inConstruct": false },
    { "name": "e1", "type": "enhancer", "score": .9, "start": 400, "end": 500, "before_start": 0, "after_end": 0, "inConstruct": false },
    { "name": "e2", "type": "enhancer", "score": .5, "start": 700, "end": 900, "before_start": 0, "after_end": 0, "inConstruct": false },
    { "name": "e3", "type": "enhacer", "score": .7, "start": 950, "end": 5000, "before_start": 0, "after_end": 0, "inConstruct": false },
  ])
  const [warnings, setWarnings] = useState(["Construct has no promoter."])

  function downloadDesign(e) {
    // TODO: build download construct functionality 
    console.log("do something and download")
  }

  function recalculateScore(e) {
    // TODO: build download construct functionality 
    console.log("do something and recalculate/update score")
  }

  useEffect(() => {
    setDesignSize(genomeSize-cargoSize)
  }, [genomeSize, cargoSize]);

  useEffect(() => {
    let newWarnings = [];
    let total_size = 0;
    let num_promoters = 0;
    regulatoryRegions.forEach(element => {
      if (element.inConstruct) {
        if (element.type === "promoter") {
          num_promoters++;
        }
        total_size = total_size + ((element.end + element.after_end) - (element.start - element.before_start));
      }
    });
    if (num_promoters > 1) {
      newWarnings.push("Construct has more than 1 promoter.")
    }
    if (num_promoters < 1) {
      newWarnings.push("Construct has no promoter.")
    }
    if (total_size > designSize) {
      newWarnings.push("Construct is: " + total_size + "bp larger then available size(" + designSize + "bp).")
    }
    setWarnings(newWarnings);
  }, [regulatoryRegions]);

  function changeConstruct(e) {
    const index = e.target.getAttribute('index');
    let newRegulatoryRegions = Array.from(regulatoryRegions);
    newRegulatoryRegions[index]["inConstruct"] = !newRegulatoryRegions[index]["inConstruct"];
    setRegulatoryRegions(newRegulatoryRegions);
  }

  function handleChange(e) {
    const type = e.target.getAttribute('objectkey');
    const index = e.target.getAttribute('index');
    const value = e.target.value;
    let newRegulatoryRegions = Array.from(regulatoryRegions);
    if (value === "") {
      newRegulatoryRegions[index][type] = 0;
    } else {
      newRegulatoryRegions[index][type] = parseInt(value);
    }
    setRegulatoryRegions(newRegulatoryRegions)
  }

  return (
    <React.Fragment>

      <div className="text-center">
        <svg width="90%" height="100">
          {/* TODO: check resizing on SVG */}
          {/* available space */}
          <Construct regulatoryRegions={regulatoryRegions} genomeSize={genomeSize} />
          <rect x="0" y="50" width={`${(designSize / genomeSize) * 100}%`} height="25"
            style={{ fill: "#defcf3" }} />
          {/* cargo size space */}
          <rect x={`${(designSize / genomeSize) * 100}%`} y="50" width="100%" height="25"
            style={{ fill: "#6bcaaf" }} />
          <text x="0" y="92" fill="black">0 bp</text>
          <text x={`${(designSize / genomeSize) * 100 - 1}%`} y="92" fill="black">{designSize} bp</text>
          <text x="96.5%" y="92" fill="black">{genomeSize} bp</text>
        </svg>
      </div>
      <hr />
      <h2>Construct Info</h2>
      <form>
        <div className="form-group row">
          <label htmlFor="inputGenome" className="col-sm-2 col-form-label">Genome Size</label>
          <div className="col-sm-5">
            <input value={genomeSize} onChange={e => setGenomeSize(e.target.value)} type="number" className="form-control" id="inputGenome" />
          </div>
        </div>
        <div className="form-group row">
          <label htmlFor="inputCargo" className="col-sm-2 col-form-label">Cargo Size</label>
          <div className="col-sm-5">
            <input value={cargoSize} onChange={e => setCargoSize(e.target.value)} type="number" className="form-control" id="inputCargo" />
          </div>
        </div>
        <div className="form-group row">
          <label htmlFor="inputConstruct" className="col-sm-2 col-form-label">Construct Size</label>
          <div className="col-sm-5">
            <input value={designSize} type="number" className="form-control" id="inputConstruct" disabled/>
          </div>
        </div>
      </form>



      <h2>Available Regions</h2>
      <div className="regionsContainer">
        <div className="regionHeaders3">Name</div>
        <div className="regionHeaders1">Type</div>
        <div className="regionHeaders3">Score</div>
        <div className="regionHeaders1">Start</div>
        <div className="regionHeaders1">End</div>
        <div className="regionHeaders2">Before Start</div>
        <div className="regionHeaders2">After End</div>
        <div className="regionHeaders3">Include</div>
        <div className="regionHeaders3">Recalculate score</div>
      </div>

      <ReactSortable list={regulatoryRegions} setList={setRegulatoryRegions}>

        {regulatoryRegions.map((item, index) => {
          return (
            <div key={index} className="regionsContainer">
              <div className="regionRows3">{item.name}</div>
              <div className="regionRows1">{item.type}</div>
              <div className="regionRows3">{item.score}</div>
              <div className="regionRows1">{item.start}</div>
              <div className="regionRows1">{item.end}</div>
              <div className="regionRows2"><input type="number" className="form-control" placeholder="bp before start"
                onChange={handleChange} objectkey="before_start" index={index} /></div>
              <div className="regionRows2"><input type="number" className="form-control" placeholder="bp after end"
                onChange={handleChange} objectkey="after_end" index={index} /></div>
              <div className="regionRows3 text-center"><input type="checkbox" onChange={changeConstruct} index={index}
                aria-label="Checkbox for following text input" /></div>
              <div className="regionRows3 text-center">
                <button type="button" onClick={recalculateScore} className="btn btn-custom2">
                  <FontAwesomeIcon icon={faSyncAlt} />

                </button>
              </div>
            </div>
          )
        })}
      </ReactSortable>

      {warnings.map((item, index) => {
        return (<div key={index} className="alert alert-warning" role="alert">
          {item}
        </div>)
      })}


      <button onClick={downloadDesign} type="button"
        className="btn btn-submit float-right">Download Design</button>

    </React.Fragment>)
}