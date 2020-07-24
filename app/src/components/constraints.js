import React, { useState } from 'react';
import ConstraintButton from './constraintButton'
import NavButtons from './navButtons';

export default function Constraints(props) {
  const [start, setStart] = useState(-1);
  const [end, setEnd] = useState(-1);
  const [n, setN] = useState(-1);
  const [select, setSelect] = useState({
    "tad": false,
    "kbregion": false,
    "gene2gene": false,
    "custom": false,
  });

  function handleClick(type){
    let newSelect = {
      "tad": false,
      "kbregion": false,
      "gene2gene": false,
      "custom": false,
    };
    newSelect[type] = true;
    setSelect(newSelect);
  }

  return (
    <React.Fragment>
      <h1 className="text-center">Set Contraints</h1>
      <hr />
      <div className="row spacing">
        <div className="col-sm-1"></div>
        <div className="col-sm-5">
          <ConstraintButton handleClick={handleClick} selectState={select["gene2gene"]} type="gene2gene" />
        </div>
        <div className="col-sm-5">
          <ConstraintButton handleClick={handleClick} selectState={select["tad"]} type="tad" />
        </div>
        <div className="col-sm-1"></div>
      </div>
      <div className="row spacing">
        <div className="col-sm-1"></div>
        <div className="col-sm-5">
          <ConstraintButton setN={setN} handleClick={handleClick} 
          selectState={select["kbregion"]} type="kbregion" />
        </div>
        <div className="col-sm-5">
          <ConstraintButton setStart={setStart} setEnd={setEnd} 
          handleClick={handleClick} selectState={select["custom"]} type="custom" />
        </div>
        <div className="col-sm-1"></div>
      </div>
      <hr />
      <NavButtons next={props.next} back={props.back} page={props.page} />
    </React.Fragment>)
}