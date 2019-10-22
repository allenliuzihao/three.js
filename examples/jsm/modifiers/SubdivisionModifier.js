/**
 *  @author Zihao Liu for Catmull-Clark subdivision and texture mapping onto the surface
 *	@author zz85 / http://twitter.com/blurspline / http://www.lab4games.net/zz85/blog
 *	@author centerionware / http://www.centerionware.com
 *
 *	Subdivision Geometry Modifier
 *		using Catmull-Clark and Loop Subdivision Scheme
 *
 *	References:
 *  	Catmull-Clark Subdivision:
 * 			http://users.cms.caltech.edu/~cs175/cs175-02/resources/CC.pdf
 * 			https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface
 *		Loop Subdivision
 *			http://graphics.stanford.edu/~mdfisher/subdivision.html
 *			http://www.holmes3d.net/graphics/subdivision/
 *			http://www.cs.rutgers.edu/~decarlo/readings/subdiv-sg00c.pdf
 *
 *	Known Issues:
 *		- currently doesn't handle "Sharp Edges" for Loop scheme
 */

import { Face3 } from '../../../src/core/Face3.js';
import { Quad } from '../../../src/core/Quad.js';
import { Geometry } from '../../../src/core/Geometry.js';
import { Vector2 } from '../../../src/math/Vector2.js';
import { Vector3 } from '../../../src/math/Vector3.js';

var SubdivisionModifier = function ( subdivisions, scheme ) {

	this.subdivisions = ( subdivisions === undefined ) ? 1 : subdivisions;

	this.scheme = scheme;

};

/**
 * constants definition
 */
SubdivisionModifier.prototype.constants = {};
SubdivisionModifier.prototype.constants.CATMULL_CLARK_SCHEME = "catmull_clark";
SubdivisionModifier.prototype.constants.LOOP_SCHEME = "loop";

// Applies the "modify" pattern
/**
 * @param {Geometry} geometry
 */
SubdivisionModifier.prototype.modify = function ( geometry ) {

	if ( geometry.isBufferGeometry ) {

		geometry = new Geometry().fromBufferGeometry( geometry );

	} else {

		geometry = geometry.clone();

	}

	geometry.mergeVertices();
	geometry.toQuadMesh();

	var repeats = this.subdivisions;

	while ( repeats -- > 0 ) {
		this.smooth( geometry );
	}

	// threejs render triangles
	geometry.toTriangleMesh();

	// compute normals
	geometry.computeFaceNormals();
	geometry.computeVertexNormals();

	return geometry;

};

( function () {

	/////////////////////////////
	// Performs one iteration of Subdivision
	SubdivisionModifier.prototype.smooth = function ( geometry ) {
		if(this.scheme === SubdivisionModifier.prototype.constants.LOOP_SCHEME) {
			Loop_Subdivision (geometry);
		} else if (this.scheme === SubdivisionModifier.prototype.constants.CATMULL_CLARK_SCHEME){
			Catmull_Clark_Subdivision (geometry);
		}
	};

	/* Catmull-Clark subdivision scheme */
	const EDGE_SPLIT = "_";

	/**
	 * Given an array of vectors, compute its average
	 * @param {Vector3[]} vectors 
	 */
	function avg(vectors){
		if(vectors.length === 0){
			return new Vector3();
		}

		let result = new Vector3();
		for (const vector of vectors) {
			result.add(vector);
		}
		return result.multiplyScalar(1/vectors.length);
	}

	/**
	 * generate edge key given vertex indices of the edge.
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 */
	function generateEdgeKey (vertexAIdx, vertexBIdx) {
		let vertexIndexA = Math.min( vertexAIdx, vertexBIdx );
		let vertexIndexB = Math.max( vertexAIdx, vertexBIdx );
		return vertexIndexA + EDGE_SPLIT + vertexIndexB;
	}

	/**
	 * generate "vertex1-vertex2" to faces set per edge
	 * 
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 * @param {Number} faceIdx 
	 * @param {Map<String, Set<Number>>} edgeToFaces 
	 */
	function generateEachEdge (vertexAIdx, vertexBIdx, faceIdx, edgeToFaces){
		let key = generateEdgeKey(vertexAIdx, vertexBIdx);
		if(!(edgeToFaces.has(key))){
			edgeToFaces.set(key, new Set([faceIdx]));
		} else {	
			edgeToFaces.get(key).add(faceIdx);
		}
	}

	/**
	 * generate "vertex1-vertex2" to faces set per face
	 * 
	 * @param {Face3[]|Quad[]} faces 
	 * @param {Map<String, Set<Number>>} edgeToFaces 
	 */
	function generateEdges(faces, edgeToFaces) {
		let face;

		for (let i = 0, il = faces.length; i < il; i++ ) {
			face = faces[i];

			if(face.constructor.name === "Face3"){
				generateEachEdge(face.a, face.b, i, edgeToFaces);
				generateEachEdge(face.b, face.c, i, edgeToFaces);
				generateEachEdge(face.c, face.a, i, edgeToFaces);
			} else if (face.constructor.name === "Quad"){
				generateEachEdge(face.a, face.b, i, edgeToFaces);
				generateEachEdge(face.b, face.c, i, edgeToFaces);
				generateEachEdge(face.c, face.d, i, edgeToFaces);
				generateEachEdge(face.d, face.a, i, edgeToFaces);
			}
		}
	}

	/**
	 * add a face information to a vertex.
	 * @param {Number} vertexIdx 
	 * @param {Number} faceIdx 
	 * @param {Map<Number, Set<Number>>} vertexToFaces 
	 */
	function addFacePerVertex (vertexIdx, faceIdx, vertexToFaces){
		let faceSet;
		if(vertexToFaces.has(vertexIdx)){
			faceSet = vertexToFaces.get(vertexIdx);
		} else {
			faceSet = new Set();
			vertexToFaces.set(vertexIdx, faceSet);
		}
		faceSet.add(faceIdx);
	}

	/**
	 * add edge information to a vertex.
	 * @param {Number} vertexAIdx 
	 * @param {Number} vertexBIdx 
	 * @param {Map<Number, Set<String>>} vertexToEdges
	 */
	function addEdgePerVertex (vertexAIdx, vertexBIdx, vertexToEdges){
		let edgeSet;
		if(vertexToEdges.has(vertexAIdx)){
			edgeSet = vertexToEdges.get(vertexAIdx);
		} else {
			edgeSet = new Set();
			vertexToEdges.set(vertexAIdx, edgeSet);
		}
		edgeSet.add(generateEdgeKey(vertexAIdx, vertexBIdx));

		if(vertexToEdges.has(vertexBIdx)){
			edgeSet = vertexToEdges.get(vertexBIdx);
		} else {
			edgeSet = new Set();
			vertexToEdges.set(vertexBIdx, edgeSet);
		}
		edgeSet.add(generateEdgeKey(vertexAIdx, vertexBIdx));
	}

	/**
	 * generate mapping of vertex -> {set {faces}, set {edges}}
	 * basically all the incident faces and edges for that vertex 
	 * 
	 * @param {Face3[]|Quad[]} faces 
	 * @param {Map<Number, Set<Number>>} vertexToFaces
	 * @param {Map<Number, Set<String>>} vertexToEdges
	 */
	function generateVertexInfo(faces, vertexToFaces, vertexToEdges){
		let face;
		for (let fi = 0, fl = faces.length; fi < fl; ++fi) {
			face = faces[fi];
			if(face.constructor.name === "Face3"){
				addFacePerVertex(face.a, fi, vertexToFaces);
				addFacePerVertex(face.b, fi, vertexToFaces);
				addFacePerVertex(face.c, fi, vertexToFaces);

				addEdgePerVertex(face.a, face.b, vertexToEdges);
				addEdgePerVertex(face.b, face.c, vertexToEdges);
				addEdgePerVertex(face.c, face.a, vertexToEdges);
			} else if (face.constructor.name === "Quad"){
				addFacePerVertex(face.a, fi, vertexToFaces);
				addFacePerVertex(face.b, fi, vertexToFaces);
				addFacePerVertex(face.c, fi, vertexToFaces);
				addFacePerVertex(face.d, fi, vertexToFaces);

				addEdgePerVertex(face.a, face.b, vertexToEdges);
				addEdgePerVertex(face.b, face.c, vertexToEdges);
				addEdgePerVertex(face.c, face.d, vertexToEdges);
				addEdgePerVertex(face.d, face.a, vertexToEdges);
			}
		}
	}

	/**
	 * calculate face point per face
	 * 
	 * @param {Vector3[]} vertices 
	 * @param {Face3[]|Quad[]} faces 
	 * @param {Map<Number, Vector3>} faceToFacePoints 
	 */
	function generateFacePoints (vertices, faces, faceToFacePoints){
		let face, vertexA, vertexB, vertexC;

		for (let index = 0, length = faces.length; index < length; index++) {
			face = faces[index];

			if (face.constructor.name === "Face3"){
				vertexA = vertices[face.a];
				vertexB = vertices[face.b];
				vertexC = vertices[face.c];
				faceToFacePoints.set(index, 
					new Vector3()
						.addVectors(vertexA, vertexB)
						.add(vertexC)
						.multiplyScalar(1/3));
			} else if (face.constructor.name === "Quad"){
				vertexA = vertices[face.a];
				vertexB = vertices[face.b];
				vertexC = vertices[face.c];
				let vertexD = vertices[face.d];

				faceToFacePoints.set(index, 
					new Vector3()
						.addVectors(vertexA, vertexB)
						.add(vertexC)
						.add(vertexD)
						.multiplyScalar(1/4));
			}
		}
	}

	/**
	 * for each edge, generate an edge point.
	 * 
	 * @param {Vector3[]} vertices 
	 * @param {Map<String, Set<Number>>} edgeToFaces
	 * @param {Map<Number, Vector3>} faceToFacePoints
	 * @param {Map<String, Vector3>} edgeToEdgePoints
	 */
	function generateEdgePoints (vertices, edgeToFaces, faceToFacePoints, edgeToEdgePoints) {
		let edgeSplits, facesArray;
		let vertexA, vertexB, facePointA, facePointB;
		let edgePoint;
	
		for (let [edge, faces] of edgeToFaces) {
			edgeSplits = edge.split(EDGE_SPLIT);
			vertexA = vertices[Number(edgeSplits[0])];
			vertexB = vertices[Number(edgeSplits[1])];

			if(faces.size !== 2){
				console.error("SubdivisionModifier generateEdgePoints error: each edge should only have two incident faces.");
				return;
			}

			facesArray = Array.from(faces);
			facePointA = faceToFacePoints.get(facesArray[0]);
			facePointB = faceToFacePoints.get(facesArray[1]);

			edgePoint = 
				new Vector3()
				.addVectors(vertexA, vertexB)
				.add(facePointA)
				.add(facePointB)
				.multiplyScalar(1/4);
			
			edgeToEdgePoints.set(edge, edgePoint);
		}
	}

	/**
	 * 
	 * @param {Vector3[]} vertices 
	 * @param {Map<Number, Set<Face3|Quad>>} vertexToFaces 
	 * @param {Map<Number, Set<String>>} vertexToEdges 
	 * @param {Map<Face3|Quad, Vector3>} faceToFacePoint
	 * @param {Vector3[]} nVertices 
	 */
	function generateNewVertexPosition (vertices, vertexToFaces, vertexToEdges, faceToFacePoint, nVertices) {
		if(nVertices.length !== vertices.length){
			return;
		}

		let vertex, avgFacePoints, avgEdgeMidPoints, n;
		let edgeSplits, vertexA, vertexB;
		let faces, edges, facePoints, edgeMidPoints;
		for (let index = 0; index < vertices.length; index++) {
			vertex = vertices[index];
			faces = vertexToFaces.get(index);
			edges = vertexToEdges.get(index);
			facePoints = [];
			edgeMidPoints = [];

			if(faces.size == 0 || edges.size == 0){
				console.error("SubdivisionModifier generateNewVertexPosition: faces and edges for vertex should not be empty.", this);
				return;
			}

			for (const face of faces) {
				if(faceToFacePoint.has(face)){
					facePoints.push(faceToFacePoint.get(face));
				} else {
					console.error("SubdivisionModifier generateNewVertexPosition: face should have a face point.", this);
					return;
				}
			}

			for (const edge of edges) {
				edgeSplits = edge.split(EDGE_SPLIT);
				vertexA = vertices[Number(edgeSplits[0])];
				vertexB = vertices[Number(edgeSplits[1])];
				edgeMidPoints.push(new Vector3().addVectors(vertexA, vertexB).multiplyScalar(1/2));
			}

			if(facePoints.length === 0 || edgeMidPoints.length === 0){
				console.error("SubdivisionModifier generateNewVertexPosition: number of face points and edge midpoints should not be empty.", this);
				return;
			}

			if(facePoints.length !== edgeMidPoints.length) {
				console.error("SubdivisionModifier generateNewVertexPosition: number of face points should match the number of edge midpoints", this);
				return;
			}

			n = facePoints.length;
			avgFacePoints = avg(facePoints);
			avgEdgeMidPoints = avg(edgeMidPoints);
			nVertices[index] = 
				new Vector3()
				.add(avgFacePoints)
				.addScaledVector(avgEdgeMidPoints, 2)
				.addScaledVector(vertex, n - 3)
				.multiplyScalar(1/n);
		}
	}

	/**
	 * @param {Number} vertexA
	 * @param {Number} vertexB
	 * @param {Number} offset 
	 * @param {Map<String, Vector3>} edgeToEdgePoint 
	 * @param {Map<Vector3, Number>} edgePointToEdgePointIndex 
	 * @param {Vector3[]} edgePoints 
	 */
	function saveEdgePoint(faceVertexA, faceVertexB, offset, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints){
		if(edgePoints === undefined || !(edgePoints instanceof Array)){
			console.error("edge points must be valid: ", edgePoints);
			return -1;
		}

		let edgePointIndex;
		let edgePoint = edgeToEdgePoint.get(generateEdgeKey(faceVertexA, faceVertexB));

		if(!edgePointToEdgePointIndex.has(edgePoint)){
			edgePointIndex = offset + edgePoints.length;
			edgePointToEdgePointIndex.set(edgePoint, edgePointIndex);
			edgePoints.push(edgePoint);
		} else {
			edgePointIndex = edgePointToEdgePointIndex.get(edgePoint);
		}

		return edgePointIndex;
	}

	/**
	 * @param {Number} faceIdx 
	 * @param {Vector2[][][]} oldUvs 
	 * @param {Vector2[][][]} newUvs 
	 * @param {boolean} quad 
	 */
	function generateNewUvs(faceIdx, oldUvs, newUvs, quad){
		if(oldUvs === undefined){
			console.warn("generateNewUvs could not find old uvs.");
			return;
		}

		if(quad){
			let uva, uvb, uvc, uvd;
			let uvfp, uve1, uve2, uve3, uve4;
			for(let l = 0, ll = oldUvs.length; l < ll; ++l){
				uva = oldUvs[l][faceIdx][0];
				uvb = oldUvs[l][faceIdx][1];
				uvc = oldUvs[l][faceIdx][2];
				uvd = oldUvs[l][faceIdx][3];
				uvfp = 
					new Vector2()
					.addVectors(uva, uvb)
					.add(uvc)
					.add(uvd)
					.multiplyScalar(1/4);
				uve1 = 
					new Vector2()
					.addVectors(uva, uvb)
					.multiplyScalar(1/2);
				uve2 = 
					new Vector2()
					.addVectors(uvb, uvc)
					.multiplyScalar(1/2);
				uve3 = 
					new Vector2()
					.addVectors(uvc, uvd)
					.multiplyScalar(1/2);
				uve4 = 
					new Vector2()
					.addVectors(uvd, uva)
					.multiplyScalar(1/2);

				newUvs[l].push([ uva, uve1, uvfp, uve4 ]);
				newUvs[l].push([ uve1, uvb, uve2, uvfp ]);
				newUvs[l].push([ uvfp, uve2, uvc, uve3 ]);
				newUvs[l].push([ uve4, uvfp, uve3, uvd ]);
			}
		} else {
			let uva, uvb, uvc;
			let uvfp, uve1, uve2, uve3;
			for(let l = 0, ll = oldUvs.length; l < ll; ++l){
				uva = oldUvs[l][faceIdx][0];
				uvb = oldUvs[l][faceIdx][1];
				uvc = oldUvs[l][faceIdx][2];
				uvfp = 
					new Vector2()
					.addVectors(uva, uvb)
					.add(uvc)
					.multiplyScalar(1/3);
				uve1 = 
					new Vector2()
					.addVectors(uva, uvb)
					.multiplyScalar(1/2);
				uve2 = 
					new Vector2()
					.addVectors(uvb, uvc)
					.multiplyScalar(1/2);
				uve3 = 
					new Vector2()
					.addVectors(uva, uvc)
					.multiplyScalar(1/2);

				newUvs[l].push([ uve1, uvb, uve2, uvfp ]);
				newUvs[l].push([ uvfp, uve2, uvc, uve3 ]);
				newUvs[l].push([ uvfp, uve3, uva, uve1 ]);
			}
		}
	}

	/**
	 * generate quad faces and new vertices from face points and edge points.
	 * 
	 * @param {Vector3[]} vertices
	 * @param {Face3[]|Quad[]} faces
	 * @param {Vector2[][][]} faceVertexUvs
	 * @param {Map<Number, Vector3>} faceToFacePoint
	 * @param {Map<String, Vector3>} edgeToEdgePoint
	 * @returns object containing vertices and quads of smoothed mesh
	 */
	function generateSmoothedMesh (vertices, faces, faceVertexUvs, faceToFacePoint, edgeToEdgePoint) {
		let quads = new Array();
		let facePoints = new Array(), edgePoints = new Array();
		let edgePointToEdgePointIndex = new Map();

		let uvs = new Array(faceVertexUvs.length);
		for(let i = 0, il = faceVertexUvs.length; i < il; ++i){
			uvs[i] = new Array();
		}

		let facePointIndex = 0;
		let face;
		for (let [faceIdx, facePoint] of faceToFacePoint) {
			face = faces[faceIdx];
			if (face.constructor.name === "Face3"){
				facePoints.push(facePoint);

				let edgePoint1Index = saveEdgePoint(face.a, face.b, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);
				let edgePoint2Index = saveEdgePoint(face.b, face.c, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);
				let edgePoint3Index = saveEdgePoint(face.c, face.a, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);

				if(edgePoint1Index === -1 || edgePoint2Index === -1 || edgePoint3Index === -1){
					return;
				}

				quads.push(
					new Quad(
						edgePoint1Index, 
						face.b, 
						edgePoint2Index, 
						facePointIndex + vertices.length,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						facePointIndex + vertices.length,
						edgePoint2Index,
						face.c, 
						edgePoint3Index,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						facePointIndex + vertices.length,
						edgePoint3Index,
						face.a, 
						edgePoint1Index, 
						undefined,
						undefined,
						face.materialIndex));

				generateNewUvs(faceIdx, faceVertexUvs, uvs, false);

				facePointIndex++;
			} else if (face.constructor.name === "Quad"){
				facePoints.push(facePoint);

				let edgePoint1Index = saveEdgePoint(face.a, face.b, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);
				let edgePoint2Index = saveEdgePoint(face.b, face.c, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);
				let edgePoint3Index = saveEdgePoint(face.c, face.d, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);
				let edgePoint4Index = saveEdgePoint(face.d, face.a, vertices.length + faceToFacePoint.size, edgeToEdgePoint, edgePointToEdgePointIndex, edgePoints);

				if(edgePoint1Index === -1 || edgePoint2Index === -1 || edgePoint3Index === -1 || edgePoint4Index === -1){
					return;
				}

				quads.push(
					new Quad(
						face.a,
						edgePoint1Index,
						facePointIndex + vertices.length, 
						edgePoint4Index,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						edgePoint1Index,
						face.b, 
						edgePoint2Index,
						facePointIndex + vertices.length,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						facePointIndex + vertices.length,
						edgePoint2Index, 
						face.c, 
						edgePoint3Index,
						undefined,
						undefined,
						face.materialIndex));
				quads.push(
					new Quad(
						edgePoint4Index,
						facePointIndex + vertices.length, 
						edgePoint3Index,
						face.d,
						undefined,
						undefined,
						face.materialIndex));

				generateNewUvs(faceIdx, faceVertexUvs, uvs, true);

				facePointIndex++;
			}
		}

		return { "vertices": vertices.concat(facePoints).concat(edgePoints), "quads": quads, "uvs": uvs};
	}

	/**
 	 * @param {Geometry} geometry
 	 */
	function Catmull_Clark_Subdivision (geometry){
		if(geometry === undefined || geometry.vertices === undefined
			|| geometry.faces == undefined || geometry.vertices.length === 0
			|| geometry.faces.length === 0){
				console.error("subdivision input geometry is not valid.", geometry);
				return;
			}

		let vertices = geometry.vertices, faces = geometry.faces, faceVertexUvs = geometry.faceVertexUvs;
		let edgeToFaces = new Map()
			,faceToFacePoint = new Map()
			,edgeToEdgePoint = new Map()
			,vertexToFaces = new Map()
			,vertexToEdges = new Map()
			,nVertices = new Array(vertices.length);

		// generate edges information; each edge information is "vertex1-vertex2" -> set {faces ids}
		generateEdges(faces, edgeToFaces);

		// generate vertex -> set {face ids} and vertex -> set {edges}
		generateVertexInfo(faces, vertexToFaces, vertexToEdges);

		// calculate the face points for all faces
		generateFacePoints(vertices, faces, faceToFacePoint);
		
		// calculate edge points
		generateEdgePoints(vertices, edgeToFaces, faceToFacePoint, edgeToEdgePoint);

		// calculate new vertex positions
		generateNewVertexPosition(vertices, vertexToFaces, vertexToEdges, faceToFacePoint, nVertices);

		// generate smoothed mesh
		let smoothed = generateSmoothedMesh(nVertices, faces, faceVertexUvs, faceToFacePoint, edgeToEdgePoint);

		// set the new geometry
		geometry.vertices = smoothed["vertices"];
		geometry.faces = smoothed["quads"];
		geometry.faceVertexUvs = smoothed["uvs"];		
	}

	/* Loop subdivision */
	
	// Some constants
	var ABC = [ 'a', 'b', 'c' ];
	
	function getEdge( a, b, map ) {

		var vertexIndexA = Math.min( a, b );
		var vertexIndexB = Math.max( a, b );

		var key = vertexIndexA + "_" + vertexIndexB;

		return map[ key ];

	}

	function processEdge( a, b, vertices, map, face, metaVertices ) {

		var vertexIndexA = Math.min( a, b );
		var vertexIndexB = Math.max( a, b );

		var key = vertexIndexA + "_" + vertexIndexB;

		var edge;

		if ( key in map ) {

			edge = map[ key ];

		} else {

			var vertexA = vertices[ vertexIndexA ];
			var vertexB = vertices[ vertexIndexB ];

			edge = {

				a: vertexA, // pointer reference
				b: vertexB,
				newEdge: null,
				// aIndex: a, // numbered reference
				// bIndex: b,
				faces: [] // pointers to face

			};

			map[ key ] = edge;

		}

		edge.faces.push( face );

		metaVertices[ a ].edges.push( edge );
		metaVertices[ b ].edges.push( edge );


	}

	/**
	 * @param {Vector3[]} vertices
	 * @param {Face3[]} faces
	 */
	function generateLookups( vertices, faces, metaVertices, edges ) {

		var i, il, face;

		for ( i = 0, il = vertices.length; i < il; i ++ ) {

			metaVertices[ i ] = { edges: [] };

		}

		for ( i = 0, il = faces.length; i < il; i ++ ) {

			face = faces[ i ];

			processEdge( face.a, face.b, vertices, edges, face, metaVertices );
			processEdge( face.b, face.c, vertices, edges, face, metaVertices );
			processEdge( face.c, face.a, vertices, edges, face, metaVertices );

		}

	}

	function newFace( newFaces, a, b, c, materialIndex ) {

		newFaces.push( new Face3( a, b, c, undefined, undefined, materialIndex ) );

	}

	function midpoint( a, b ) {

		return ( Math.abs( b - a ) / 2 ) + Math.min( a, b );

	}

	function newUv( newUvs, a, b, c ) {

		newUvs.push( [ a.clone(), b.clone(), c.clone() ] );

	}

	/**
 	 * @param {Geometry} geometry
 	 */
	function Loop_Subdivision (geometry){
		var tmp = new Vector3();

		var oldVertices, oldFaces, oldUvs;
		var newVertices, newFaces, newUVs = [];

		var n, i, il, j, k;
		var metaVertices, sourceEdges;

		// new stuff.
		var newEdgeVertices, newSourceVertices;
		
		oldVertices = geometry.vertices; // { x, y, z}
		oldFaces = geometry.faces; // { a: oldVertex1, b: oldVertex2, c: oldVertex3 }
		oldUvs = geometry.faceVertexUvs[ 0 ];

		var hasUvs = oldUvs !== undefined && oldUvs.length > 0;

		/******************************************************
		 *
		 * Step 0: Preprocess Geometry to Generate edges Lookup
		 *
		 *******************************************************/

		metaVertices = new Array( oldVertices.length );
		sourceEdges = {}; // Edge => { oldVertex1, oldVertex2, faces[]  }

		generateLookups( oldVertices, oldFaces, metaVertices, sourceEdges );


		/******************************************************
		 *
		 *	Step 1.
		 *	For each edge, create a new Edge Vertex,
		 *	then position it.
		 *
		 *******************************************************/

		newEdgeVertices = [];
		var other, currentEdge, newEdge, face;
		var edgeVertexWeight, adjacentVertexWeight, connectedFaces;

		for ( i in sourceEdges ) {

			currentEdge = sourceEdges[ i ];
			newEdge = new Vector3();

			edgeVertexWeight = 3 / 8;
			adjacentVertexWeight = 1 / 8;

			connectedFaces = currentEdge.faces.length;

			// check how many linked faces. 2 should be correct.
			if ( connectedFaces != 2 ) {

				// if length is not 2, handle condition
				edgeVertexWeight = 0.5;
				adjacentVertexWeight = 0;

				if ( connectedFaces != 1 ) {

					// console.warn( 'Subdivision Modifier: Number of connected faces != 2, is: ', connectedFaces, currentEdge );

				}

			}

			newEdge.addVectors( currentEdge.a, currentEdge.b ).multiplyScalar( edgeVertexWeight );

			tmp.set( 0, 0, 0 );

			for ( j = 0; j < connectedFaces; j ++ ) {

				face = currentEdge.faces[ j ];

				for ( k = 0; k < 3; k ++ ) {

					other = oldVertices[ face[ ABC[ k ] ] ];
					if ( other !== currentEdge.a && other !== currentEdge.b ) break;

				}

				tmp.add( other );

			}

			tmp.multiplyScalar( adjacentVertexWeight );
			newEdge.add( tmp );

			currentEdge.newEdge = newEdgeVertices.length;
			newEdgeVertices.push( newEdge );

			// console.log(currentEdge, newEdge);

		}

		/******************************************************
		 *
		 *	Step 2.
		 *	Reposition each source vertices.
		 *
		 *******************************************************/

		var beta, sourceVertexWeight, connectingVertexWeight;
		var connectingEdge, connectingEdges, oldVertex, newSourceVertex;
		newSourceVertices = [];

		for ( i = 0, il = oldVertices.length; i < il; i ++ ) {

			oldVertex = oldVertices[ i ];

			// find all connecting edges (using lookupTable)
			connectingEdges = metaVertices[ i ].edges;
			n = connectingEdges.length;

			if ( n == 3 ) {

				beta = 3 / 16;

			} else if ( n > 3 ) {

				beta = 3 / ( 8 * n ); // Warren's modified formula

			}

			// Loop's original beta formula
			// beta = 1 / n * ( 5/8 - Math.pow( 3/8 + 1/4 * Math.cos( 2 * Math. PI / n ), 2) );

			sourceVertexWeight = 1 - n * beta;
			connectingVertexWeight = beta;

			if ( n <= 2 ) {

				// crease and boundary rules
				// console.warn('crease and boundary rules');

				if ( n == 2 ) {

					// console.warn( '2 connecting edges', connectingEdges );
					sourceVertexWeight = 3 / 4;
					connectingVertexWeight = 1 / 8;

					// sourceVertexWeight = 1;
					// connectingVertexWeight = 0;

				} else if ( n == 1 ) {

					// console.warn( 'only 1 connecting edge' );

				} else if ( n == 0 ) {

					// console.warn( '0 connecting edges' );

				}

			}

			newSourceVertex = oldVertex.clone().multiplyScalar( sourceVertexWeight );

			tmp.set( 0, 0, 0 );

			for ( j = 0; j < n; j ++ ) {

				connectingEdge = connectingEdges[ j ];
				other = connectingEdge.a !== oldVertex ? connectingEdge.a : connectingEdge.b;
				tmp.add( other );

			}

			tmp.multiplyScalar( connectingVertexWeight );
			newSourceVertex.add( tmp );

			newSourceVertices.push( newSourceVertex );

		}


		/******************************************************
		 *
		 *	Step 3.
		 *	Generate Faces between source vertices
		 *	and edge vertices.
		 *
		 *******************************************************/

		newVertices = newSourceVertices.concat( newEdgeVertices );
		var sl = newSourceVertices.length, edge1, edge2, edge3;
		newFaces = [];

		var uv, x0, x1, x2;
		var x3 = new Vector2();
		var x4 = new Vector2();
		var x5 = new Vector2();

		for ( i = 0, il = oldFaces.length; i < il; i ++ ) {

			face = oldFaces[ i ];

			// find the 3 new edges vertex of each old face

			edge1 = getEdge( face.a, face.b, sourceEdges ).newEdge + sl;
			edge2 = getEdge( face.b, face.c, sourceEdges ).newEdge + sl;
			edge3 = getEdge( face.c, face.a, sourceEdges ).newEdge + sl;

			// create 4 faces.

			newFace( newFaces, edge1, edge2, edge3, face.materialIndex );
			newFace( newFaces, face.a, edge1, edge3, face.materialIndex );
			newFace( newFaces, face.b, edge2, edge1, face.materialIndex );
			newFace( newFaces, face.c, edge3, edge2, face.materialIndex );

			// create 4 new uv's

			if ( hasUvs ) {

				uv = oldUvs[ i ];

				x0 = uv[ 0 ];
				x1 = uv[ 1 ];
				x2 = uv[ 2 ];

				x3.set( midpoint( x0.x, x1.x ), midpoint( x0.y, x1.y ) );
				x4.set( midpoint( x1.x, x2.x ), midpoint( x1.y, x2.y ) );
				x5.set( midpoint( x0.x, x2.x ), midpoint( x0.y, x2.y ) );

				newUv( newUVs, x3, x4, x5 );
				newUv( newUVs, x0, x3, x5 );

				newUv( newUVs, x1, x4, x3 );
				newUv( newUVs, x2, x5, x4 );

			}

		}

		// Overwrite old arrays
		geometry.vertices = newVertices;
		geometry.faces = newFaces;
		if ( hasUvs ) geometry.faceVertexUvs[ 0 ] = newUVs;

		// console.log('done');
	}

} )();

export { SubdivisionModifier };
