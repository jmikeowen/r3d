/*
 *
 *  	v2d.c
 *
 *  	See v2d.h for usage.
 *
 *  	Devon Powell
 *  	15 October 2015
 *
 *		Copyright (c) 2015, The Board of Trustees of the Leland Stanford Junior University, 
 *		through SLAC National Accelerator Laboratory (subject to receipt of any required approvals 
 *		from the U.S. Dept. of Energy). All rights reserved. Redistribution and use in source and 
 *		binary forms, with or without modification, are permitted provided that the following 
 *		conditions are met: 
 *		(1) Redistributions of source code must retain the above copyright notice, 
 *		this list of conditions and the following disclaimer. 
 *		(2) Redistributions in binary form must reproduce the above copyright notice, 
 *		this list of conditions and the following disclaimer in the documentation and/or other 
 *		materials provided with the distribution. 
 *		(3) Neither the name of the Leland Stanford Junior University, SLAC National Accelerator 
 *		Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to endorse 
 *		or promote products derived from this software without specific prior written permission. 
 *		
 *		THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
 *		OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 *		MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
 *		COPYRIGHT OWNER, THE UNITED STATES GOVERNMENT, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
 *		INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 *		LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
 *		BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
 *		STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
 *		USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *		
 *		You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to 
 *		the features, functionality or performance of the source code ("Enhancements") to anyone; 
 *		however, if you choose to make your Enhancements available either publicly, or directly to 
 *		SLAC National Accelerator Laboratory, without imposing a separate written license agreement 
 *		for such Enhancements, then you hereby grant the following license: a non-exclusive, 
 *		royalty-free perpetual license to install, use, modify, prepare derivative works, 
 *		incorporate into other computer software, distribute, and sublicense such Enhancements or 
 *		derivative works thereof, in binary and source code form.
 *
 */

#include "v2d.h"
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define wav(va, wa, vb, wb, vr) {			\
	vr.x = (wa*va.x + wb*vb.x)/(wa + wb);	\
	vr.y = (wa*va.y + wb*vb.y)/(wa + wb);	\
}

// TODO: make this a generic "split" routine that just takes a plane.
void r2d_split(r2d_poly* inpoly, r2d_poly** outpolys, r2d_real coord, r2d_int ax);

void r2d_rasterize(r2d_poly* poly, r2d_dvec2 ibox[2], r2d_real* dest_grid, r2d_rvec2 d, r2d_int polyorder) {

	r2d_int i, m, spax, dmax, nstack, siz;
	r2d_int nmom = R2D_NUM_MOMENTS(polyorder);
	r2d_real moments[nmom];
	r2d_poly* children[2];
	r2d_dvec2 gridsz;

	// return if any parameters are bad 
	for(i = 0; i < 2; ++i) gridsz.ij[i] = ibox[1].ij[i]-ibox[0].ij[i];	
	if(!poly || poly->nverts <= 0 || !dest_grid || 
			gridsz.i <= 0 || gridsz.j <= 0) return;
	
	// explicit stack-based implementation
	// stack size should never overflow in this implementation, 
	// even for large input grids (up to ~512^2) 
        m = (r2d_int)(ceil(log2(gridsz.i))+ceil(log2(gridsz.j))+1);
	struct {
		r2d_poly poly;
		r2d_dvec2 ibox[2];
	} stack[m];
        for (i = 0; i < m; ++i) stack[i].poly = r2d_init_empty_poly();

	// push the original polyhedron onto the stack
	// and recurse until child polyhedra occupy single rasters
	nstack = 0;
	r2d_copy_poly(&stack[nstack].poly, poly);
	memcpy(stack[nstack].ibox, ibox, 2*sizeof(r2d_dvec2));
	nstack++;
	while(nstack > 0) {

		// pop the stack
		// if the leaf is empty, skip it
		--nstack;
		if(stack[nstack].poly.nverts <= 0) continue;
		
		// find the longest axis along which to split 
		dmax = 0;
		spax = 0;
		for(i = 0; i < 2; ++i) {
			siz = stack[nstack].ibox[1].ij[i]-stack[nstack].ibox[0].ij[i];
			if(siz > dmax) {
				dmax = siz; 
				spax = i;
			}	
		}

		// if all three axes are only one raster long, reduce the single raster to the dest grid
#define gind(ii, jj, mm) (nmom*((ii-ibox[0].i)*gridsz.j+(jj-ibox[0].j))+mm)
		if(dmax == 1) {
			r2d_reduce(&stack[nstack].poly, moments, polyorder);
			// TODO: cell shifting for accuracy
			for(m = 0; m < nmom; ++m)
				dest_grid[gind(stack[nstack].ibox[0].i, stack[nstack].ibox[0].j, m)] += moments[m];
			continue;
		}

		// split the poly and push children to the stack
		children[0] = &stack[nstack].poly;
		children[1] = &stack[nstack+1].poly;
		r2d_split(&stack[nstack].poly, children, d.xy[spax]*(stack[nstack].ibox[0].ij[spax]+dmax/2), spax);
		memcpy(stack[nstack+1].ibox, stack[nstack].ibox, 2*sizeof(r2d_dvec2));
		stack[nstack].ibox[1].ij[spax] -= dmax-dmax/2; 
		stack[nstack+1].ibox[0].ij[spax] += dmax/2;
		nstack += 2;
	}
}

void r2d_split(r2d_poly* inpoly, r2d_poly** outpolys, r2d_real coord, r2d_int ax) {

	// direct access to vertex buffer
	if(inpoly->nverts <= 0) return;
	r2d_int* nverts = &inpoly->nverts;
	r2d_int v, np, onv, vcur, vnext, vstart, nright, cside;
	r2d_rvec2 newpos;
	r2d_int* side = malloc((*nverts)*sizeof(r2d_int));
	r2d_real sdists[(*nverts)];
	r2d_vertex* verts_new;
	r2d_int* side_new;

	printf("Starting:--------------------------------------------------------------------------------\n");
	printf("inpoly %d\n", inpoly->nverts);
	r2d_print(inpoly);
	printf("outpolys[0] %d\n", outpolys[0]->nverts);
	r2d_print(outpolys[0]);
	printf("outpolys[1] %d\n", outpolys[1]->nverts);
	r2d_print(outpolys[1]);

	// calculate signed distances to the clip plane
	nright = 0;
	for (v = 0; v < *nverts; ++v) side[v] = 0;
	for(v = 0; v < *nverts; ++v) {
			//sdists[v] = splane.d + r2d_dot(inpoly->verts[v].pos, splane.n);
		sdists[v] = coord - inpoly->verts[v].pos.xy[ax];
		sdists[v] *= -1;
		if(sdists[v] < 0.0) {
			side[v] = 1;
			nright++;
		}
	}

	// return if the poly lies entirely on one side of it 
	if(nright == 0) {
                r2d_copy_poly(outpolys[0], inpoly);
		outpolys[1]->nverts = 0;
		return;
	}
	if(nright == *nverts) {
                r2d_copy_poly(outpolys[1], inpoly);
		outpolys[0]->nverts = 0;
		return;
	}

	// check all edges and insert new vertices on the bisected edges 
	onv = inpoly->nverts;
	for(vcur = 0; vcur < onv; ++vcur) {
		if(side[vcur]) continue;
		for(np = 0; np < 2; ++np) {
			vnext = inpoly->verts[vcur].pnbrs[np];
			if(!side[vnext]) continue;
			wav(inpoly->verts[vcur].pos, -sdists[vnext],
				inpoly->verts[vnext].pos, sdists[vcur],
				newpos);
                        verts_new = realloc(inpoly->verts, (*nverts + 2)*sizeof(r2d_vertex));
			side_new = realloc(side, (*nverts + 2)*sizeof(r2d_int));
			if (verts_new == NULL || side_new == NULL) {
                          printf("r2d_split ERROR: unable to increase vertex buffer size.\n");
                          abort();
                        }
                        inpoly->verts = verts_new;
			side = side_new;
			side[*nverts] = 0;
			inpoly->verts[*nverts].pos = newpos;
			inpoly->verts[vcur].pnbrs[np] = *nverts;
			inpoly->verts[*nverts].pnbrs[np] = -1;
			inpoly->verts[*nverts].pnbrs[1-np] = vcur;
			(*nverts)++;
			side[*nverts] = 1;
			inpoly->verts[*nverts].pos = newpos;
			inpoly->verts[*nverts].pnbrs[1-np] = -1;
			inpoly->verts[*nverts].pnbrs[np] = vnext;
			inpoly->verts[vnext].pnbrs[1-np] = *nverts;
			(*nverts)++;
		}
	}
	printf("--------------------------------------------------------------------------------\n");
	r2d_print(&inpoly);

	abort();

	// for each new vert, search around the poly for its new neighbors
	// and doubly-link everything
	for(vstart = onv; vstart < *nverts; ++vstart) {
		if(inpoly->verts[vstart].pnbrs[1] >= 0) continue;
		vcur = inpoly->verts[vstart].pnbrs[0];
		do {
			vcur = inpoly->verts[vcur].pnbrs[0]; 
		} while(vcur < onv);
		inpoly->verts[vstart].pnbrs[1] = vcur;
		inpoly->verts[vcur].pnbrs[0] = vstart;
	}
	printf("--------------------------------------------------------------------------------\n");
	r2d_print(&inpoly);

	// copy and compress vertices into their new buffers
	// reusing side[] for reindexing
	onv = *nverts;
	outpolys[0]->nverts = 0;
	outpolys[1]->nverts = 0;
	outpolys[0]->verts = (r2d_vertex*) malloc(onv*sizeof(r2d_vertex));
	outpolys[1]->verts = (r2d_vertex*) malloc(onv*sizeof(r2d_vertex));
	if (outpolys[0]->verts == NULL || outpolys[1]->verts == NULL) {
		printf("r2d_split ERROR: unable to allocate outpoly vertex buffers.\n");
		abort();
	}
	for(v = 0; v < onv; ++v) {
		cside = side[v];
		if (cside != 0 && cside != 1) {
			printf("NOPE: %d %d\n", v, cside);
		}
		outpolys[cside]->verts[outpolys[cside]->nverts] = inpoly->verts[v];
		side[v] = (outpolys[cside]->nverts)++;
	}

	for(v = 0; v < outpolys[0]->nverts; ++v) 
		for(np = 0; np < 2; ++np)
			outpolys[0]->verts[v].pnbrs[np] = side[outpolys[0]->verts[v].pnbrs[np]];
	for(v = 0; v < outpolys[1]->nverts; ++v) 
		for(np = 0; np < 2; ++np)
			outpolys[1]->verts[v].pnbrs[np] = side[outpolys[1]->verts[v].pnbrs[np]];
	free(side);
}

void r2d_get_ibox(r2d_poly* poly, r2d_dvec2 ibox[2], r2d_rvec2 d) {
	r2d_int i, v;
	r2d_rvec2 rbox[2];
	for(i = 0; i < 2; ++i) {
		rbox[0].xy[i] = 1.0e30;
		rbox[1].xy[i] = -1.0e30;
	}
	for(v = 0; v < poly->nverts; ++v) {
		for(i = 0; i < 2; ++i) {
			if(poly->verts[v].pos.xy[i] < rbox[0].xy[i]) rbox[0].xy[i] = poly->verts[v].pos.xy[i];
			if(poly->verts[v].pos.xy[i] > rbox[1].xy[i]) rbox[1].xy[i] = poly->verts[v].pos.xy[i];
		}
	}
	for(i = 0; i < 2; ++i) {
		ibox[0].ij[i] = floor(rbox[0].xy[i]/d.xy[i]);
		ibox[1].ij[i] = ceil(rbox[1].xy[i]/d.xy[i]);
	}
}

void r2d_clamp_ibox(r2d_poly* poly, r2d_dvec2 ibox[2], r2d_dvec2 clampbox[2], r2d_rvec2 d) {
	r2d_int i, nboxclip;
	r2d_plane boxfaces[4];
	nboxclip = 0;
	memset(boxfaces, 0, sizeof(boxfaces));
	for(i = 0; i < 2; ++i) {
		if(ibox[1].ij[i] <= clampbox[0].ij[i] || ibox[0].ij[i] >= clampbox[1].ij[i]) {
			memset(ibox, 0, sizeof(ibox));
			poly->nverts = 0;
			return;
		}
		if(ibox[0].ij[i] < clampbox[0].ij[i]) {
			ibox[0].ij[i] = clampbox[0].ij[i];
			boxfaces[nboxclip].d = -clampbox[0].ij[i]*d.xy[i];
			boxfaces[nboxclip].n.xy[i] = 1.0;
			nboxclip++;
		}
		if(ibox[1].ij[i] > clampbox[1].ij[i]) {
			ibox[1].ij[i] = clampbox[1].ij[i];
			boxfaces[nboxclip].d = clampbox[1].ij[i]*d.xy[i];
			boxfaces[nboxclip].n.xy[i] = -1.0;
			nboxclip++;
		}	
	}
	if(nboxclip) r2d_clip(poly, boxfaces, nboxclip);
}
