#pragma once

enum EHungarianMatchingType {
	HUNGARIAN_MATCH_MIN = 0,
	HUNGARIAN_MATCH_MAX = 1,
};

/*
 * findMatching
 * implementation of the Hungarian matching algorithm
 * referenced from: http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
 */
bool HungarianMatching(const MatrixReal& m, MatrixReal& result, EHungarianMatchingType type, int maxsteps = 10000);
