import React from 'react';
import NavButtons from './navButtons';
import IGV from './igv'

export default function Genes(props) {

  return (
    <React.Fragment>
      <h1 className="text-center">Gene Selection</h1>
      <hr />
      <div className="row">
        <div className="col-sm-12">
          <h3>Search</h3>
          <div class="input-group">
            <input type="text" class="form-control" placeholder="gene symbol" />
            <div class="input-group-append" id="button-addon4">
              <button class="btn btn-outline-custom" type="button">Search for Genes</button>
              <button class="btn btn-outline-custom" type="button">Suggest Genes</button>
            </div>
          </div>
          <hr />
          <h3>Gene List</h3>
          <div className="form-group">
            <select multiple className="form-control" id="exampleFormControlSelect2">
            </select>
          </div>
        </div>
      </div>
      <IGV genome={'hg19'} start={50000} end={500000} chrom={'chr1'} />
        <NavButtons next={props.next} back={props.back} page={props.page} />
    </React.Fragment>)
}