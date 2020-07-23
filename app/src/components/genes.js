import React from 'react';


export default function Genes(props) {

  return (
    <React.Fragment>
      <h1 className="text-center">Gene Selection</h1>
      <hr />
      <div className="row">
        <div className="col-sm-12">
          <h3>Sample List</h3>
          <div className="form-group">
            <select multiple className="form-control" id="exampleFormControlSelect2">
            </select>
          </div>
        </div>
      </div>
    </React.Fragment>)
}