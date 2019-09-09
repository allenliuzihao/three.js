import { Vector3 } from '../math/Vector3';
import { Color } from '../math/Color';

export interface Event {
	type: string;
	target?: any;
	[attachment: string]: any;
}

/**
 * Quad face.
 *
 * @source Quad.js
 */
export class Quad {

	/**
	 * @param a Vertex A index.
	 * @param b Vertex B index.
	 * @param c Vertex C index.
	 * @param d Vertex D index.
	 * @param normal Face normal or array of vertex normals.
	 * @param color Face color or array of vertex colors.
	 * @param materialIndex Material index.
	 */
	constructor(
		a: number,
		b: number,
		c: number,
		d: number,
		normal?: Vector3,
		color?: Color,
		materialIndex?: number
	);
	constructor(
		a: number,
		b: number,
		c: number,
		d: number,
		normal?: Vector3,
		vertexColors?: Color[],
		materialIndex?: number
	);
	constructor(
		a: number,
		b: number,
		c: number,
		d: number,
		vertexNormals?: Vector3[],
		color?: Color,
		materialIndex?: number
	);
	constructor(
		a: number,
		b: number,
		c: number,
		d: number,
		vertexNormals?: Vector3[],
		vertexColors?: Color[],
		materialIndex?: number
	);

	/**
	 * Vertex A index.
	 */
	a: number;

	/**
	 * Vertex B index.
	 */
	b: number;

	/**
	 * Vertex C index.
	 */
	c: number;

	/**
	 * Vertex D index.
	 */
	d: number;

	/**
	 * Face normal.
	 */
	normal: Vector3;

	/**
	 * Array of 4 vertex normals.
	 */
	vertexNormals: Vector3[];

	/**
	 * Face color.
	 */
	color: Color;

	/**
	 * Array of 4 vertex normals.
	 */
	vertexColors: Color[];

	/**
	 * Material index (points to {@link Geometry.materials}).
	 */
	materialIndex: number;
	
	toTriangleFaces(): Face3[];
	toQuadMesh(): Quad;

	clone(): this;

	copy( source: Quad ): this;

}
