import React from 'react';

function RegionList(props) {
    // props
    // regions
    // type
    // enzymes
    // tfs
    const regs = props.regions.filter(reg =>
        reg.type===props.type
        ).map((r, index) =>
            <tr key={props.type + index}>
                <td>{r.chrom}</td>
                <td>{r.strand}</td>
                <td>{r.start}</td>
                <td>{r.end}</td>
                <td>{r.score.toFixed(3)}</td>
                <td>
                <button type="button" className="btn btn-secondary" onClick={e => props.selectRegion(index, props.type)}>Select region
                </button>
                </td>
            </tr>
        );

    return (
        <tbody>
            {regs}
        </tbody>
    )
}

export default RegionList;
