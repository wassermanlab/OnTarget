import React from 'react';
import Nav from '../components/nav'

export default function About () {
    return (
        <React.Fragment>
            <Nav />
            <div className="container navbar-fix">
            <h1>About</h1>
            <p>
            All of our cells have the same DNA sequences (genes) that “code” for 
            the body’s entire diversity of cell types, from blood to brain, liver, 
            muscle and skin. The differences between cell types are determined 
            by which genes are active, and at what times, in each cell. 
            Our chromosomes include thousands of scattered DNA sequences 
            called regulatory regions (or On/Off switches) that control 
            when, where, and to what degree, genes are turned on/off. 
            These “On/Off switches” function differently in each cell type. 
            Changes (i.e. mutations) in the DNA sequence of an On/Off switch 
            can lead to abnormal gene activity, which, in turn, can cause disease.
            <br/><br/>
            Within the OnTarget project, we have created new computer programs to analyze vast amounts of data and:
            <br/>
            - identify these On/Off switches in each cell type
            <br/>
            - predict the effect of mutations in the DNA sequences of these On/Off switches
            - design DNA sequences to turn a gene on/off for therapy.
            <br/><br/>
            The OnTarget project has built upon the success of databases and 
            software created in our Canadian lab. Potential benefits to 
            Canadians include improved methods to treat patients with genetic 
            disorders and advanced training of a new generation of highly-skilled 
            data scientists.
            </p>
            <h1>Creators</h1>
            <p>
                Oriol Fornes<br/>
                Tamar Av-Shalom<br/>
            </p>
            <h1>Video Tutorial</h1>
            <p></p>
            <h1>Manuscript</h1>
            <p></p>
            <h1>Code and Data</h1>
            <p></p>



            </div>
        </React.Fragment>)
}