import React, { useState } from 'react';
import './app.css';
import Results from '../components/results'
import Design from '../components/design'
import Nav from '../components/nav'
import { useParams } from "react-router-dom";

export default function Res() {
  const [tab, setTab] = useState("design")
  const id = useParams().id;

  function rendered_component() {
    switch (tab) {
      case "landscape":
        return (<Results id={id} />);
      case "design":
        return (<Design id={id} />);

      default:
        return (<React.Fragment></React.Fragment>);
    }
  }

  return (
    <React.Fragment>
      <Nav />
      <div className="navbar-fix container-fluid">
        <div className='result-id'><b>Result ID:</b> {id}</div>


        <ul className="nav nav-tabs">
          <li className="nav-item">
            <div className={`${tab === "landscape" ? "active" : ""} nav-link`} onClick={() => setTab("landscape")}>
              Regulatory Landscape</div>
          </li>
          <li className="nav-item">
            <div className={`${tab === "design" ? "active" : ""} nav-link`} onClick={() => setTab("design")}>
              Design Module</div>
          </li>
        </ul>
        {rendered_component()}
      </div>
    </React.Fragment>)
}