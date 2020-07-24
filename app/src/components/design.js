import React, { useState, forwardRef } from 'react';
import NavButtons from './navButtons';
import IGV from './igv'
import { ReactSortable } from "react-sortablejs";


const CustomComponent = forwardRef((props, ref) => {
  return (

    <table className="table">
      <thead>
        <tr>
          <th scope="col">Name</th>
          <th scope="col">Type</th>
          <th scope="col">Core Region Score</th>
          <th scope="col">Start</th>
          <th scope="col">End</th>
          <th scope="col">bp before start</th>
          <th scope="col">bp after end</th>
        </tr>
      </thead>
      <tbody>
        {props.children}
      </tbody>
    </table>
  )
});


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

  function checkWarnings(regs) {
    // more than 1 promoter, construct too long, regions overlap eachother
    // todo add region overlap eachother check
    let newWarnings = [];
    let total_size = 0;
    let num_promoters = 0;
    regs.forEach(element => {
      console.log(element)
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
  }

  function changeConstruct(e) {
    const index = e.target.getAttribute('index');
    let newRegulatoryRegions = Array.from(regulatoryRegions);
    newRegulatoryRegions[index]["inConstruct"] = !newRegulatoryRegions[index]["inConstruct"];
    checkWarnings(newRegulatoryRegions);
    setRegulatoryRegions(newRegulatoryRegions);
  }

  function handleChange(e) {
    const type = e.target.getAttribute('key');
    const index = e.target.getAttribute('index');
    const value = e.target.value;
    let newRegulatoryRegions = Array.from(regulatoryRegions);
    newRegulatoryRegions[index][type] = value;
    setRegulatoryRegions(newRegulatoryRegions)
  }

  return (
    <React.Fragment>
      <h1 className="text-center">Design</h1>
      <hr />

      <svg width="100%" height="50">
        <g> {/* available space */}
          <rect x="0" y="0" width={`${(designSize / genomeSize) * 100}%`} height="100%"
            style={{ fill: "white", stroke: "black", strokeWidth: 1, }} />
        </g>
        <g> {/* cargo */}
          <rect x={`${(designSize / genomeSize) * 100}%`} y="0"
            width={`${(cargoSize / genomeSize) * 100}%`} height="100%"
            style={{ fill: "black", stroke: "black", strokeWidth: 1, }} />
          <line x1={`${(designSize / genomeSize) * 100}%`}
            y1="0" x2={`${(designSize / genomeSize) * 100}%`} y2="100%"
            style={{ stroke: "green", strokeWidth: 3 }} />
        </g>
      </svg>

      <hr />

      <h2>Available Regions</h2>
      <div>
        <div class="input-group mb-3">
          <div class="input-group-prepend">
            <span class="input-group-text">name</span>
            <span class="input-group-text">type</span>
            <span class="input-group-text">score</span>
            <span class="input-group-text">start</span>
            <span class="input-group-text">end</span>
          </div>
          <div class="input-group-append">
            <span class="input-group-text">bp before start</span>
            <span class="input-group-text">bp after end</span>
            <span class="input-group-text">include in construct</span>
          </div>

        </div>
      </div>

      <ReactSortable list={regulatoryRegions} setList={setRegulatoryRegions}>

        {regulatoryRegions.map((item, index) => {
          return (
            <div>
              <div class="input-group mb-3">
                {/* stuff before */}
                <div class="input-group-prepend">
                  <span class="input-group-text">{item.name}</span>
                  <span class="input-group-text">{item.type}</span>
                  <span class="input-group-text">{item.score}</span>
                  <span class="input-group-text">{item.start}</span>
                  <span class="input-group-text">{item.end}</span>
                </div>
                {/* text inputs */}
                <input type="number" class="form-control" placeholder="bp before start"
                  onChange={handleChange} key="before_start" index={index} />
                <input type="number" class="form-control" placeholder="bp after end"
                  onChange={handleChange} key="after_end" index={index} />
                {/* checkmark at end */}
                <div class="input-group-append">
                  <div class="input-group-text">
                    <input type="checkbox" onChange={changeConstruct} index={index}
                      aria-label="Checkbox for following text input" />
                  </div>
                </div>
              </div>
            </div>
          )
        })}
      </ReactSortable>

      {warnings.map(item => {
        return (<div class="alert alert-warning" role="alert">
          {item}
        </div>)
      })}


      <button onClick={downloadDesign} type="button"
        className="btn btn-submit">Download Design</button>

      <hr />
      <NavButtons next={props.next} back={props.back} page={props.page} />
    </React.Fragment>)
}