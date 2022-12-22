import * as React from "react";
import { Routes, Route } from "react-router-dom";
import Layout from "./layout";
import Home from "./home";
import Design from "./design";
import NoMatch from "./noMatch";
import GetRegions from "./getRegions";


export default function App() {
  return (
    <div>
      <Routes>
        <Route path="/" element={<Layout />}>
          <Route index element={<Home />} />
          <Route path="design" element={<Design />} />
          <Route path="getregions" element={<GetRegions />} />
          <Route path="*" element={<NoMatch />} />
        </Route>
      </Routes>
    </div>
  );
}
