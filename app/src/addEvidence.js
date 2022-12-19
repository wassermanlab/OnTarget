import React from 'react';

function AddEvidence(props) {
    if (props.page !== 2) {
        return null;
    }


    return (
        <div>
            <form>
                <h2>Upload Evidence</h2>
                <div className="form-group">
                    <input type="file" multiple />
                </div>
                <div className="form-group">
                    <button className="btn btn-secondary" type="submit">Upload</button>
                </div>
            </form>

        </div>
    );





}

export default AddEvidence;
