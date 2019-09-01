import {
  BufferGeometry,
  Geometry
} from '../../../src/Three';

export class SubdivisionModifier {
  constructor(subdivisions?: number, scheme?: string);
  subdivisions: number;
  scheme: string;

  modify(geometry: BufferGeometry | Geometry): Geometry;
  smooth(geometry: Geometry): void;
}
