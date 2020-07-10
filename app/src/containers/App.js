import React, { useState } from 'react';
import './app.css';


export default function App() {
  const [page, setPageState] = useState(2)
  const [form, setForm] = useState({
    samples: [],
    gene: "",
    searchWindow: "",
  })
  const pages ={
    1: {width: "25%", step:"SAMPLES"},
    2: {width: "50%", step:"GENE"},
    3: {width: "75%", step:"CONSTRAINS"},
    4: {width: "100%", step:"RESULTS"},
  }

  return (
    <React.Fragment>
      <div className="container">
        <div className="row">
          <div className="col-sm">
            <div class="progress" style={{height:"30px"}}>
              <div class="progress-bar" role="progressbar" style={{width:pages[page]["width"]}}>{pages[page]["step"]}</div>
            </div>
          </div>
        </div>
      </div>
    </React.Fragment>)
}