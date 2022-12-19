import React from 'react';

function CreateMiniPromoter(props) {
    if (props.page !== 4) {
        return null;
    }

    const regs = props.regulatoryRegions.map((r, index) =>
        <tr key={index}>
            <td>{r.chrom}</td>
            <td>{r.strand}</td>
            <td>{r.start}</td>
            <td>{r.end}</td>
            <td>{r.type}</td>
            <td>{r.score}</td>
            <td>
                <input name="selected" checked={r.selected} onChange={e => props.handleRegCheckBox(index, e)} type="checkbox" />
            </td>
        </tr>
    );



    const miniPromoters = props.miniPromoterCart.map((miniPromoter, index) =>
        <tr key={index}>
            <td>{miniPromoter.regulatoryRegionsList.toString()}</td>
            <td>{miniPromoter.totalScore}</td>
            <td>
                <button type="button" className="btn btn-primary ontarget-button" onClick={e => props.removeMiniPromoter(index)}>Remove
                </button>
            </td>
        </tr>
    );

    return (
        <React.Fragment>
            {/* regulatoryRegions */}
            <div className="row">
                <h2>Regulatory Regions</h2>
                <table className="table">
                    <thead>
                        <tr>
                            <th scope="col">Chromosome</th>
                            <th scope="col">Strand</th>
                            <th scope="col">Start</th>
                            <th scope="col">End</th>
                            <th scope="col">Type</th>
                            <th scope="col">Score</th>
                            <th scope="col">Select Region</th>
                        </tr>
                    </thead>
                    <tbody>
                        {regs}
                    </tbody>
                </table>
            </div>
            {/* add mini promoters */}
            <div className="row row-margin">
                <button type="button" className="btn btn-primary ontarget-button" onClick={props.addMiniPromoter}>Add Mini Promoter</button>
            </div>
            {/* mini promoter card */}
            <div className="row row-margin">
                <h2>Mini Promoter Cart</h2>
                <table className="table">
                    <thead>
                        <tr>
                            <th scope="col">Regulatory Regions</th>
                            <th scope="col">Total Score</th>
                            <th scope="col"></th>
                        </tr>
                    </thead>
                    <tbody>
                        {miniPromoters}
                    </tbody>
                </table>
            </div>
            {/* download mini promoters */}
            <div className="row">
                <button type="button" className="btn btn-primary ontarget-button" onClick={props.downloadMiniPromoters}>Download Mini Promoter</button>
            </div>
        </React.Fragment>
    )
}

export default CreateMiniPromoter;
