import React, { useState } from 'react';
import './app.css';
import Genes from '../components/genes'
import Constraints from '../components/constraints'
import Samples from '../components/samples'
import Results from '../components/results'


export default function App() {
  const [page, setPageState] = useState(1)
  const [form, setForm] = useState({
    samples: [],
    gene: "",
    searchWindow: "",
  })
  function back(){
    setPageState(page-1)
  }

  function next(){
    setPageState(page+1)
  }

  function rendered_component() {
    switch (page) {
      case 1:
        return (<Samples />);
      case 2:
        return (<Genes />);
      case 3:
        return (<Constraints />);
      case 4:
        return (<Results />);
      default:
        return (<React.Fragment></React.Fragment>);
    }
  }

  return (
    <React.Fragment>
      <div className="container-fluid">
      {rendered_component()}
      <div className="btn-group float-right" role="group" aria-label="Basic example">
        {page!==1? <button type="button" onClick={back} className="btn btn-nav">Back</button> : <React.Fragment></React.Fragment>}
        {page!==4? <button type="button" onClick={next} className="btn btn-nav">Next</button> : <React.Fragment></React.Fragment>}
      </div>
      </div>
    </React.Fragment>)
}