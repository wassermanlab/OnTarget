import React from 'react';

function SelectedRegionList(props) {
    const selectedRegs = props.selectedRegions.map((r, index) =>
        <tr key={index}>
            <td>{r.type}</td>
            <td>{r.chrom}</td>
            <td>{r.strand}</td>
            <td>{r.start}</td>
            <td>{r.end}</td>
            <td>{r.score.toFixed(3)}</td>
            <td>
                <button type="button" className="btn btn-secondary" onClick={e => props.removeSelectedRegion(index)}>Select region
                </button>
            </td>
        </tr>
    );

    return (
        <tbody>
            {selectedRegs}
        </tbody>
    )
}

export default SelectedRegionList;