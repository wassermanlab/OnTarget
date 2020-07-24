import React, { useState } from 'react';
import './app.css';
import Genes from '../components/genes'
import Constraints from '../components/constraints'
import Samples from '../components/samples'
import Results from '../components/results'
import Nav from '../components/nav'

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
        return (<Samples next={next} back={back} page={page} />);
      case 2:
        return (<Genes next={next} back={back} page={page} />);
      case 3:
        return (<Constraints next={next} back={back} page={page} />);
      case 4:
        return (<Results next={next} back={back} page={page} />);
      default:
        return (<React.Fragment></React.Fragment>);
    }
  }

  return (
    <React.Fragment>
      <Nav />
      <div className="navbar-fix container-fluid">
      {rendered_component()}
        </div>
    </React.Fragment>)
}