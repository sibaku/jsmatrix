import {
    MatF32 as mf32, // shorter name, since we only use float values
    VecF32 as vf32, // shorter name, since we only use float values
} from '../../jsmatrix.js';

import * as m from '../../jsmatrix.js';


// near and far value for projection
var near = 1;
var far = 20;
// number of random points
var numPoints = 2000;

// canvas and context
var canvas;
var ctx;

// Draw a single point
function drawPoint(p) {
    // Throw away points behind the camera
    if (p.at(2) < 0) {
        return;
    }
    // We linearize the projective depth to get a simple depth effect
    var z_ndc = 2.0 * p.at(2) - 1.0;
    var z_eye = 2.0 * near * far / (far + near - z_ndc * (far - near));

    // Set color and alpha so it will be less visible further away, thus helping the missing depth test
    var factor = Math.max(0.0, Math.min(1.0, z_eye / (far - near)));
    var col = Math.max(0, Math.min(255, Math.floor(factor * factor * 255)));
    var fillStyle = "rgba(" + col + "," + col + "," + col + "," + Math.max(0.1, 1 - factor) + ")";
    ctx.fillStyle = fillStyle;
    // draw the point
    ctx.beginPath();
    ctx.arc(p.at(0), p.at(1), 4, 0, Math.PI * 2);
    ctx.fill();
}
// Draw the points
function drawPoints(ctx, mvp, viewportMat, points) {
    // Compute the full transform from local to screen
    var VP = m.mult(viewportMat, mvp);
    var p2 = m.mult(VP, points);
    // Homogenize points after projection
    m.homogenize(p2);

    // Sim
    ctx.save();

    m.colwise(p2, drawPoint);
    ctx.restore();

}

// Draw a line from one point to another
function drawLine(ctx, mvp, viewportMat, p0, p1) {
    // The same as for the full point array
    var VP = m.mult(viewportMat, mvp);
    p0 = m.mult(VP, p0);
    m.homogenize(p0);
    p1 = m.mult(VP, p1);
    m.homogenize(p1);

    // Draw a line in red from the projected first point to the second
    ctx.save();
    ctx.strokeStyle = "red";
    ctx.lineWidth = 4;
    ctx.beginPath();
    ctx.moveTo(p0.at(0), p0.at(1));
    ctx.lineTo(p1.at(0), p1.at(1));
    ctx.stroke();
    ctx.restore();
}

window.onload = function () {

    canvas = document.getElementById("canvas");
    ctx = canvas.getContext("2d");

    // Define a view matrix
    // As with OpenGL we use y as the up vector
    var V = m.lookAt(vf32.from([8, 5, 0]), vf32.from([0, 0, 0]), vf32.from([0, 1, 0]));
    // Define a perspective projection
    var P = m.perspective(m.deg2rad(90), canvas.width / canvas.height, near, far);

    // No model matrix here, so model matrix is just Projection * View
    var MVP = m.mult(P, V);
    // Viewport
    // Last parameter flips the y axis, since the transform by default assumes the origin in the lower left corner.
    // HTML canvas starts at the top left.
    var vp = m.viewport(0, 0, canvas.width, canvas.height, true);

    // Create an uninitialized matrix to store the points in
    // Typed arrays will still have zeros at the beginning, but general array may have not
    // This avoids the explicit initialization
    // Points are stored per column with a homogeneous coordinate
    var points = mf32.uninitialized(4, numPoints);

    // We want to create random points on a plane with some uncertainty
    // These points are then used to estimate a normal

    // Views for the x,y and z rows
    var rx = m.row(points, 0);
    var ry = m.row(points, 1);
    var rz = m.row(points, 2);

    // Fill each coordinate with their own range of random values
    m.insert(rx, mf32.rand(1, points.cols()));
    // Scale random [0,1] to [0,10]
    m.scale(rx, 10, rx);
    // Center random values to [-5,5]
    m.addScalar(rx, -5);

    m.insert(rz, mf32.rand(1, points.cols()));
    // Scale random [0,1] to [0,6]
    m.scale(rz, 6, rz);
    // Center random values to [-3,3]
    m.addScalar(rz, -5);

    // The plane is in the xz-plane
    // The y component will not have as much random range, just a bit to simulate a bit of noise
    m.insert(ry, mf32.rand(1, points.cols()));
    // Scale random [0,1] to [0,0.4]
    m.scale(ry, 0.4, ry);
    // Center random values to [-0.2,0.2]
    m.addScalar(rz, -0.2);

    // Set the homogeneous coordinate 1 for all points
    m.fill(m.row(points, 3), 1.0);


    // We apply a random translation and rotation
    // To compose the transforms, we use the homogeneous axis angle version

    // We use a random axis, which could theoretically be zero, but we ignore that here
    // Don't rotate too much
    var R4 = m.axisAngle4(vf32.rand(3), Math.random() * 2 * Math.PI / 10);
    // Add a random vector in [-0.5,0.5]^3
    var T = m.translation(m.addScalar(vf32.rand(3), -0.5));

    var transform = m.mult(T, R4);

    // Transform points
    // We can also apply the transform in place
    m.colwise(points, col => m.insert(col, m.mult(transform, col)));

    // To compute the normal, we first center all points, that is, subtract the center
    // For each axis (each row), we compute the sum and divide it by the number of points
    // This can be accomplished with a row reduction
    var mean = m.rowreduce(points, row => m.sum(row) / row.cols());

    // We copy the points to change them for the normal computation
    var pcenterd = m.copy(points);

    // Subtract the mean from
    m.colwise(pcenterd, c => {
        // Each column is a vector and we take the first 3 components
        let c3 = m.subvec(c, 0, 3);
        // Subtract the mean (ignoring the homogeneous coordinate) from the point in the column and store the result back in the matrix
        m.sub(c3, m.subvec(mean, 0, 3), c3);
    });



    // Given zero-centered points stored as matrix columns, we can estimate a normal
    // The normal is found as the eigenvector for the least eigenvalue of the outer product points*points^T
    // These are just the singular values of the point matrix (No details here)
    // We compute the singular value decomposition of the x,y,z coordinates. For that we call the computeSVD with a block view
    var svd = m.computeSVD(m.block(pcenterd, 0, 0, 3, pcenterd.cols()));

    // The normal is the left singular vector for the least singular value
    // The singular values are sorted from largest to smallest, so we get the last column
    var N = m.col(svd.U, 2);

    // As the direction is not well specified, we flip it according to the y-axis
    // As long as we are not perpendicular to that, this will mostly give us normals pointing in roughly the same direction
    if (m.dot(N, m.Y) < 0.0) {
        m.neg(N, N);
    }

    // For drawing purposes, we will create a second point to draw a line from the mean point to
    var endpoint = m.copy(mean);
    var endpoint3 = m.subvec(endpoint, 0, endpoint.rows() - 1);
    // Add the scaled normal to the mean to get the final endpoint
    m.add(endpoint3, m.scale(N, 2), endpoint3);

    // Draw the points
    drawPoints(ctx, MVP, vp, points);

    // Draw the normal located at the mean
    drawLine(ctx, MVP, vp, mean, endpoint);
};