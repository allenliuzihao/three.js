import { Color } from '../math/Color.js';
import { Vector3 } from '../math/Vector3.js';
import { Face3 } from './Face3.js';

/**
 * @author Zihao Liu
 */

function Quad( a, b, c, d, normal, color, materialIndex ) {

	this.a = a;
	this.b = b;
	this.c = c;
	this.d = d;

	this.normal = ( normal && normal.isVector3 ) ? normal : new Vector3();
	this.vertexNormals = Array.isArray( normal ) ? normal : [];

	this.color = ( color && color.isColor ) ? color : new Color();
	this.vertexColors = Array.isArray( color ) ? color : [];

	this.materialIndex = materialIndex !== undefined ? materialIndex : 0;

}

Object.assign( Quad.prototype, {

	clone: function () {

		return new this.constructor().copy( this );

	},

	copy: function ( source ) {

		this.a = source.a;
		this.b = source.b;
		this.c = source.c;
		this.d = source.d;

		this.normal.copy( source.normal );
		this.color.copy( source.color );

		this.materialIndex = source.materialIndex;

		for ( var i = 0, il = source.vertexNormals.length; i < il; i ++ ) {

			this.vertexNormals[ i ] = source.vertexNormals[ i ].clone();

		}

		for ( i = 0, il = source.vertexColors.length; i < il; i ++ ) {

			this.vertexColors[ i ] = source.vertexColors[ i ].clone();

		}

		return this;

	},

	/**
	 * counter-clockwise from a to d
	 */
	toTriangleFaces: function (){
		let face1 = new Face3(this.a, this.b, this.d, this.normal, this.color, this.materialIndex);
		if (this.vertexColors !== undefined && this.vertexColors.length === 4){
			face1.vertexColors.push(this.vertexColors[0]);
			face1.vertexColors.push(this.vertexColors[1]);
			face1.vertexColors.push(this.vertexColors[3]);
		}
		if (this.vertexNormals !== undefined && this.vertexNormals.length === 4){
			face1.vertexNormals.push(this.vertexNormals[0]);
			face1.vertexNormals.push(this.vertexNormals[1]);
			face1.vertexNormals.push(this.vertexNormals[3]);
		}

		let face2 = new Face3(this.b, this.c, this.d, this.normal, this.color, this.materialIndex);
		if (this.vertexColors !== undefined && this.vertexColors.length === 4){
			face2.vertexColors.push(this.vertexColors[1]);
			face2.vertexColors.push(this.vertexColors[2]);
			face2.vertexColors.push(this.vertexColors[3]);
		}
		if (this.vertexNormals !== undefined && this.vertexNormals.length === 4){
			face2.vertexNormals.push(this.vertexNormals[1]);
			face2.vertexNormals.push(this.vertexNormals[2]);
			face2.vertexNormals.push(this.vertexNormals[3]);
		}
		
		return [face1, face2];
	}

} );


export { Quad };
