import * as React from "react";
import { Routes, Route } from "react-router-dom";
import Layout from "./layout";
import Home from "./home";
import Design from "./design";
import NoMatch from "./noMatch";
import Example from "./example";
import ADORA2Example from "./adora2example";
import "../node_modules/bootstrap/dist/css/bootstrap.min.css";
import "../node_modules/bootstrap/dist/js/bootstrap.bundle.min.js";
//import GetRegions from "./getRegions";


export default function App() {
  return (
    <div>
      <Routes>
        <Route path="/" element={<Layout />}>
          <Route index element={<Home />} />
          <Route path="design" element={<Design />} />
          <Route path="pitx3example" element={<Example />} />
          <Route path="adora2example" element={<ADORA2Example />} />
          <Route path="*" element={<NoMatch />} />
        </Route>
      </Routes>
    </div>
  );
}
