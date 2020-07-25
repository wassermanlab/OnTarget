import React, { useState, useEffect } from 'react';
import NavButtons from './navButtons';
import IGV from './igv'
import { ReactSortable } from "react-sortablejs";


function Construct(props) {
  // regulatoryRegions, genomeSize
  const [c, setC] = useState([])


  useEffect(() => {
    let offset = 0;
    let newC = []
    props.regulatoryRegions.forEach(e=>{
      if(e.inConstruct) {
        const lenPercentage = (((e.end + e.after_end) - (e.start - e.before_start))/props.genomeSize)*100;
        const y2 = offset + lenPercentage;
        newC.push({"type": e.type, "name": e.name, "y1": offset, "y2": y2})
        offset = y2;
      }
    })
    setC(newC);
  }, [props.regulatoryRegions, props.genomeSize]);
  let offset = 0;
  return (
    <g>
    {c.map((item, index) => {
      return (<line x1={`${item.y1}%`} y1="20" x2={`${item.y2}%`} 
      y2="20" style={{ stroke: (item.type==="promoter"? "orange": "blue"), strokeWidth: 10 }} />)
    })}
    </g>)
}


export default function Design(props) {
  const [genomeSize, setGenomeSize] = useState(4500);
  const [cargoSize, setCardoSize] = useState(2500);
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
    console.log("do something and download")
  }

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
    newRegulatoryRegions[index][type] = parseInt(value);
    setRegulatoryRegions(newRegulatoryRegions)
  }

  return (
    <React.Fragment>
      <h1 className="text-center">Design</h1>
      <hr />

      <svg width="100%" height="100">

        {/* available space */}
        <Construct regulatoryRegions={regulatoryRegions} genomeSize={genomeSize} />
        <line x1="0" y1="90" x2={`${(designSize / genomeSize) * 100}%`} y2="90"
          style={{ stroke: "green", strokeWidth: 10 }} />
        {/* available space */}
        <line x1={`${(designSize / genomeSize) * 100}%`} y1="90" x2="100%" y2="90"
          style={{ stroke: "red", strokeWidth: 10 }} />
      </svg>

      <hr />

      <h2>Available Regions</h2>
      <div>
        <div className="input-group mb-3">
          <div className="input-group-prepend">
            <span className="input-group-text">name</span>
            <span className="input-group-text">type</span>
            <span className="input-group-text">score</span>
            <span className="input-group-text">start</span>
            <span className="input-group-text">end</span>
          </div>
          <div className="input-group-append">
            <span className="input-group-text">bp before start</span>
            <span className="input-group-text">bp after end</span>
            <span className="input-group-text">include in construct</span>
          </div>

        </div>
      </div>

      <ReactSortable list={regulatoryRegions} setList={setRegulatoryRegions}>

        {regulatoryRegions.map((item, index) => {
          return (
            <div key={index}>
              <div className="input-group mb-3">
                {/* stuff before */}
                <div className="input-group-prepend">
                  <span className="input-group-text">{item.name}</span>
                  <span className="input-group-text">{item.type}</span>
                  <span className="input-group-text">{item.score}</span>
                  <span className="input-group-text">{item.start}</span>
                  <span className="input-group-text">{item.end}</span>
                </div>
                {/* text inputs */}
                <input type="number" className="form-control" placeholder="bp before start"
                  onChange={handleChange} objectkey="before_start" index={index} />
                <input type="number" className="form-control" placeholder="bp after end"
                  onChange={handleChange} objectkey="after_end" index={index} />
                {/* checkmark at end */}
                <div className="input-group-append">
                  <div className="input-group-text">
                    <input type="checkbox" onChange={changeConstruct} index={index}
                      aria-label="Checkbox for following text input" />
                  </div>
                </div>
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
        className="btn btn-submit">Download Design</button>

      <hr />
      <NavButtons next={props.next} back={props.back} page={props.page} />
    </React.Fragment>)
}