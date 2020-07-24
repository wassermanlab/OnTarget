import React from 'react';
import './constraintButton.css';

export default function ConstraintButton(props) {
  const info = {
    "kbregion": { title: "n Kilobase Region", img: "/kb_around_button.svg" },
    "custom": { title: "Custom", img: null },
    "gene2gene": { title: "Gene to Gene", img: "/gene2gene_button.svg" },
    "tad": { title: "TAD Boundry", img: "/tad_button.svg" },
  }


  if (props.type === "custom") {
    return (
      <div onClick={()=> {props.handleClick(props.type)}} className={`text-center c-button${props.selectState ? "-selected" : ""}`}>
        <span className={`${props.selectState ? "dot-selected" : "dot"}`} ></span>
        <h3>{info[props.type]["title"]}</h3>
        <img style={{ width: "100%" }} src={info[props.type]["img"]} alt="" />
        <div style={{ paddingTop: "18%" }} className="input-group">
          <div className="input-group-prepend">
            {/* TODO: change chr to whatever chromosome selected gene is on */}
            <span className="input-group-text">chrX</span>
          </div>
          <input onChange={(e)=>(props.setStart(e.target.value))} type="text" placeholder="start" className="form-control" />
          <input onChange={(e)=>(props.setEnd(e.target.value))} type="text" placeholder="end" className="form-control" />
        </div>
      </div>)
  }
  else if (props.type === "kbregion") {
    return (
      <div onClick={()=> {props.handleClick(props.type)}} className={`text-center c-button${props.selectState ? "-selected" : ""}`}>
        <span className={`${props.selectState ? "dot-selected" : "dot"}`} ></span>
        <h3>{info[props.type]["title"]}</h3>
        <img style={{ width: "100%" }} src={info[props.type]["img"]} alt="" />
        {/* TODO: change max holding */}
        <input onChange={(e)=>(props.setN(e.target.value))} style={{ float:"right", width: "50px" }} className="form-control" type="number" min="1" max="5" placeholder="n" />
      </div>)
  } else {
    return (
      <div onClick={()=> {props.handleClick(props.type)}} className={`text-center c-button${props.selectState ? "-selected" : ""}`} >
        <span className={`${props.selectState ? "dot-selected" : "dot"}`} ></span>
        <h3>{info[props.type]["title"]}</h3>
        <img style={{ width: "100%" }} src={info[props.type]["img"]} alt="" />
      </div>)
  }
}