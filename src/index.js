// src/index.js

import * as Annotator from './C6-Annotator.js';
import * as Gene from './C6-Gene.js';
import * as Oligos from './C6-Oligos.js';
import * as PCA from './C6-PCA.js';
import * as Seq from './C6-Seq.js';
import * as Sim from './C6-Sim.js';
import * as Utils from './C6-Utils.js';

const C6 = {
  ...Annotator,
  ...Gene,
  ...Oligos,
  ...PCA,
  ...Seq,
  ...Sim,
  ...Utils
};

C6.VERSION = '1.0.11';

export default C6;
