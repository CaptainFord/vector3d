using UnityEngine;
using System.Text;
using System;

public struct Boundsd
{
	public static void scratchGrounds(){
//		DoNegativeExtentsTest();
//		DoPrecisionTest();
//		FindSmallestPossibleIncrements();

	}

	static void DoNegativeExtentsTest ()
	{
		//	This test proved that Unity's bounds implementation allowed negative extents (and negative sizes).
		//	Which means that the minimum is allowed to be greater than the maximum.
		//	I have since verified that this has the expected effect of making ALL of the subsequent calculations
		//	completely screwy.
		//	But a simple naive implementation isn't enough to account for all the discrepancies. I need to
		//	test what happens with NaN and infinity values.

		Vector3 v3 = new Vector3();
		
		Bounds b = new Bounds();
		b.center = Vector3.zero;
		b.extents = new Vector3(1,1,1);
		Bounds before = b;
		Bounds after = b;
		after.max = new Vector3(-2,1,1);
		Debug.Log ("Before: " + AllVectorsToString(before) + 
		           "\nAfter: " + AllVectorsToString(after));
	}

	static void DoPrecisionTest ()
	{
		Bounds extremeTest = new Bounds(new Vector3(1e-8f,1e12f,2.4e-9f),new Vector3(1e8f,1e-12f,-1f));
		Bounds extremeTest2 = new Bounds();
		extremeTest2.SetMinMax(new Vector3(2.38e-8f,3.14e-12f,6.969e-24f),new Vector3(2.38e8f,3.14e12f,6.969e24f));
		Debug.Log ("Center Precision Test:" + AllVectorsToString(extremeTest)
		           + "\nMin&Max Precision Test:" + AllVectorsToString(extremeTest2));
		// VERIFIED - Unity's bounds is stored as center & extents, as the precision on min was lost but the precision on center was not.
		// 		This is slightly disturbing to me ... it's odd, isn't it? Sure, it has greater precision, 
		//	but a box defined by the min & max would correspond exactly to actual Vector3 points. So it's
		//	... unaligned.
		//	How would unity transform a bounding box from local to world coordinates anyway?
		//	Well, anyways, I don't care, since this is how I'd implement it once I get to using decimals anyway. 
		//	Still seems odd.
		//	It's ALSO odd that unity doesn't check to see if the extents are negative. But that makes sense too.
		//		Oh hey ... in retrospect, I could have inferred that bounds is stored as center&extents from the 
		//	description of SetMinMax(). If it were stored as min&max, then SetMinMax() wouldn't be any faster.
		//		I DO wanna know why unity's bounds exposes only properties, and why center and extents aren't just
		//	fields if it doesn't do any value checking. Maybe there's a reason for that.
	}

	static void FindSmallestPossibleIncrements(){
		//	This just helped with implementing BumpUp and BumpDown
		int factorFp, factorFn, factorDp, factorDn;
		float valueFp, valueFn;
		double valueDp, valueDn;

		FindSmallestPossibleIncrement(out factorFp, out valueFp, true);
		FindSmallestPossibleIncrement(out factorFn, out valueFn, false);
		FindSmallestPossibleIncrement(out factorDp, out valueDp, true);
		FindSmallestPossibleIncrement(out factorDn, out valueDn, false);

		Debug.Log ("Fp: 2^-" + factorFp + " = " + valueFp + " -> " + (1f + valueFp).ToString("G30")
		           + "\nFn: 2^-" + factorFn + " = " + valueFn + " -> " + (1f + valueFn).ToString("G30")
		           + "\nDp: 2^-" + factorDp + " = " + valueDp + " -> " + (1f + valueDp).ToString("G30")
		           + "\nDn: 2^-" + factorDn + " = " + valueDn + " -> " + (1f + valueDn).ToString("G30")
		           );
	}

	static void FindSmallestPossibleIncrement(out int factor, out float value, bool add){
		factor = 0;
		float fvalue = add ? 1f : -1f;
		while(true){
			++factor;
			float fresult = (0.5f * fvalue) + 1f;
			fresult = fresult - 1f;
			if(fresult == 0){
				value = fvalue;
				return;
			} else {
				fvalue *= 0.5f;
			}
		}
	}
	static void FindSmallestPossibleIncrement(out int factor, out double value, bool add){
		factor = 0;
		double fvalue = add ? 1d : -1d;
		while(true){
			++factor;
			double fresult = (0.5d * fvalue) + 1d;
			fresult = fresult - 1d;
			if(fresult == 0){
				value = fvalue;
				return;
			} else {
				fvalue *= 0.5d;
			}
		}
	}

	public static string AllVectorsToString (Bounds b)
	{
		string f = "G7";
		return "c:" + b.center.ToString(f) + ";e:" + b.extents.ToString(f) + ";m:" + b.min.ToString(f) + ";M:" + b.max.ToString(f) + ";s:" + b.size.ToString(f);
	}
	public static string AllVectorsToString (Boundsd b)
	{
		string f = "G7";
		return "c:" + b.center.ToString(f) + ";e:" + b.extents.ToString(f) + ";m:" + b.min.ToString(f) + ";M:" + b.max.ToString(f) + ";s:" + b.size.ToString(f);
	}
	
	public const string defaultFormat = "{center}\u00B1{extents}";
	public const string allValuesFormat = "{center}\u00B1{extents}|{min}-{max}";
	public const string defaultNumberFormat = "G5";


	public static string ToString(Bounds b){
		return ToString (b, defaultFormat, defaultNumberFormat);
	}
	
	public static string ToString(Bounds b, string format){
		if(format.Contains("{min}") || format.Contains("{max}")
		   || format.Contains("{center}") || format.Contains("{extents}") 
		   || format.Contains("{size}")) {
			return ToString (b, format, defaultNumberFormat);
		} else {
			return ToString (b, defaultFormat, format);
		}
	}

	public static string ToString(Bounds b, string format, string numberFormat){
		StringBuilder result = new StringBuilder();
		int lastEnd = 0;
		int nextStart = format.IndexOf('{');
		while(nextStart >= 0) {
			int nextEnd = format.IndexOf('}', nextStart + 1);
			if(nextEnd < 0){
				break;
			}
			string tagName = format.Substring(nextStart + 1, nextEnd - nextStart - 1);
			string valueToString = ValueToString(b, tagName, numberFormat);
			if(valueToString != null) {
				result.Append(format, lastEnd, nextStart - lastEnd);
				result.Append(valueToString);
				lastEnd = nextEnd + 1;
			}
			nextStart = format.IndexOf('{', nextEnd + 1);
		}
		result.Append(format, lastEnd, format.Length - lastEnd);
		return result.ToString();
	}

	static string ValueToString (Bounds b, string tagName, string numberFormat)
	{
		switch(tagName){
		case "center":
			return b.center.ToString(numberFormat);
		case "extents":
			return b.extents.ToString(numberFormat);
		case "size":
			return b.size.ToString(numberFormat);
		case "min":
			return b.min.ToString(numberFormat);
		case "max":
			return b.max.ToString(numberFormat);
		default:
			return null;
		}
	}

	public Vector3d center;
	internal Vector3d m_extents;
	//	I think most of the operations would be faster if min and max weren't calculated...
	//	But using center+extents allows greater precision ... sort of ...
	//	Conversion from center+extents to 
//	internal Vector3d m_min, m_max;

//	public Vector3d center {
//		get {
//			return (m_min + m_max) * 0.5d;
//		}
//		set {
//
//		}
//	}
	public Vector3d extents {
		get {
			return m_extents;
		}
		set {
			m_extents = value;
		}
	}
	
	public Vector3d min {
		get {
			return center - m_extents;
		}
		set {
			Vector3d currentMax = this.max;
			this.center = (value + currentMax) * 0.5d;
			this.size = (currentMax - value);
		}
	}
	public Vector3d max {
		get {
			return center + m_extents;
		}
		set {
			Vector3d currentMin = this.min;
			this.center = (value + currentMin) * 0.5d;
			this.size = (value - currentMin);
		}
	}
	public Vector3d size {
		get {
			return m_extents * 2;
		}
		set {
			m_extents = value * 0.5;
		}
	}

	public Boundsd(Bounds b){
		this.center = (Vector3d)b.center;
		this.m_extents = (Vector3d)b.extents;
	}

	public override string ToString(){
		return ToString (defaultFormat, defaultNumberFormat);
	}
	
	public string ToString(string format){
		if(format.Contains("{min}") || format.Contains("{max}")
		   || format.Contains("{center}") || format.Contains("{extents}") 
		   || format.Contains("{size}")) {
			return ToString (format, defaultNumberFormat);
		} else {
			return ToString (defaultFormat, format);
		}
	}

	//	The format string uses named {}-delimeted arguments
	//	The names are the same as the field/property names
	public string ToString(string format, string numberFormat){
		StringBuilder result = new StringBuilder();
		int lastEnd = 0;
		int nextStart = format.IndexOf('{');
		while(nextStart >= 0) {
			int nextEnd = format.IndexOf('}', nextStart + 1);
			if(nextEnd < 0){
				break;
			}
			string tagName = format.Substring(nextStart + 1, nextEnd - nextStart - 1);
			string valueToString = ValueToString(tagName, numberFormat);
			if(valueToString != null) {
				result.Append(format, lastEnd, nextStart - lastEnd);
				result.Append(valueToString);
				lastEnd = nextEnd + 1;
			}
			nextStart = format.IndexOf('{', nextEnd + 1);
		}
		result.Append(format, lastEnd, format.Length - lastEnd);
		return result.ToString();
	}

	string ValueToString (string tagName, string numberFormat)
	{
		switch(tagName){
		case "center":
			return center.ToString(numberFormat);
		case "extents":
			return extents.ToString(numberFormat);
		case "size":
			return size.ToString(numberFormat);
		case "min":
			return min.ToString(numberFormat);
		case "max":
			return max.ToString(numberFormat);
		default:
			return null;
		}
	}

	public bool Contains (Vector3d point)
	{
		Vector3d diff = point - center;
		return Mathd.Abs (diff.x) <= m_extents.x 
				&& Mathd.Abs (diff.y) <= m_extents.y 
				&& Mathd.Abs (diff.z) <= m_extents.z;
	}

	public void Encapsulate (Vector3d point)
	{
//		EncapsulateAxis (ref center.x, ref m_extents.x, point.x);
//		EncapsulateAxis (ref center.y, ref m_extents.y, point.y);
//		EncapsulateAxis (ref center.z, ref m_extents.z, point.z);
		EncapsulateNaively(point);
	}

	void EncapsulateNaively(Vector3d point){
		EncapsulateAxisNaively (ref center.x, ref m_extents.x, point.x);
		EncapsulateAxisNaively (ref center.y, ref m_extents.y, point.y);
		EncapsulateAxisNaively (ref center.z, ref m_extents.z, point.z);
//		for(int i=0; i<3; ++i){
//			double diff = point[i] - center[i];
//			if(diff < 0d){
//				if(-diff > extents[i]){
//					double oldMax = center[i] + extents[i];
//					center[i] = (oldMax + point[i]) * 0.5d;
//					extents[i] = (oldMax - point[i]) * 0.5d;
//				}
//			} else {
//				if(diff > extents[i]){
//					double oldMin = center[i] - extents[i];
//					center[i] = (oldMin + point[i]) * 0.5d;
//					extents[i] = (point[i] - oldMin) * 0.5d;
//				}
//			}
//		}
	}
	void EncapsulateAxisNaively (ref double centerX, ref double extentsX, double pointX)
	{
		//	I think I figured it out. Any time this would result in a negative extent
		//	it instead adopts the point's value and sets the extent to zero.
		//	But the way I'm doing this, there's no point where that check makes sense ... except ...
		//	... no, there's just no place it makes sense.
		//	Wait ... unless ... unless it's defined somewhere that a bounds with negative extents
		//	sort of represents a bounding box which hasn't had any points loaded in yet.
		//	...well...if that's the case ... then why only check the NEW extents?

		//	There is nothing about this that I don't find bizarre and inexplicable.

		//	Also it seems like it might do this with NaN results, too.

		//	diff = NaN
		//	NaN > extentsX ... yes or no?
		//	if it fails, then ... of course. If I reversed the comparison, then ... no wait ... it wouldn't matter.
		//	This *incidentally* ignores NaN values. I don't actually see a reason to change that.

		double diff = pointX - centerX;
		if(Mathd.Abs(diff) > extentsX){
			double newExtentsX = (Mathd.Abs(diff) + extentsX) * 0.5d;
			if(newExtentsX < 0d){
				centerX = pointX;
				extentsX = 0d;
			} else {
				centerX += (diff < 0d ? (diff + extentsX) : (diff - extentsX)) * 0.5d;
				extentsX = newExtentsX;
			}
	    }
	}
	void EncapsulateAxisNaively0 (ref double centerX, ref double extentsX, double pointX)
	{
		//	I think I figured it out. Any time this would result in a negative extent
		//	it instead adopts the point's value and sets the extent to zero.
		//	But the way I'm doing this, there's no point where that check makes sense ... except ...
		//	... no, there's just no place it makes sense.

		double diff = pointX - centerX;
		if(diff < 0d){
			if(-diff > extentsX){
				double delta = (-diff - extentsX) * 0.5d;
				centerX -= delta;
				extentsX += delta;
			}
		} else {
			if(diff > extentsX){
				double delta = (diff - extentsX) * 0.5d;
				centerX += delta;
				extentsX += delta;
			}
		}
	}
	void EncapsulateNaively0(Vector3d point){
		//	This version produces zero extents MORE often
		Vector3d newMax = this.max;
		Vector3d newMin = this.min;

		for(int i=0; i<3; ++i){
			if(point[i] > newMax[i]){
				newMax[i] = point[i];
			}
			if(point[i] < newMin[i]){
				newMin[i] = point[i];
			}
		}
		SetMinMax(newMin, newMax);
	}

	void EncapsulateAxis (ref double centerX, ref double extentsX, double pointX)
	{
		//	Do it naively, assume extents is positive. Or don't.
		//	No, I can't do that. If I do that, it will always fail.
		double diffX = pointX - centerX;
//		if(Mathd.Abs (diffX) > Mathd.Abs(extentsX)) {
			if(extentsX > 0d) {
				if(diffX > 0d) {
					if(diffX > extentsX) {
						//	Maximum changes, minimum doesn't.
						double nexMaxX = pointX;
						double newMinX = centerX - extentsX;
						centerX = (newMinX + nexMaxX) * 0.5d;
						extentsX = (nexMaxX - newMinX) * 0.5d;
					}
				} else {
					if(-diffX > extentsX) {
						//	Then we need to expand the minimum to encompass the point
						//	The maximum stays the same
						double nexMaxX = centerX + extentsX;
						double newMinX = pointX;
						centerX = (newMinX + nexMaxX) * 0.5d;
						extentsX = (nexMaxX - newMinX) * 0.5d;
					}
				}
			} else {
				if(diffX > 0d) {
					if(diffX > -extentsX) {
						//	Maximum changes, minimum doesn't.
						double nexMaxX = pointX;
						double newMinX = centerX + extentsX;
						centerX = (newMinX + nexMaxX) * 0.5d;
						extentsX = (nexMaxX - newMinX) * -0.5d;
					}
				} else {
					if(-diffX > -extentsX) {
						//	Then we need to expand the minimum to encompass the point
						//	The maximum stays the same
						double nexMaxX = centerX - extentsX;
						double newMinX = pointX;
						centerX = (newMinX + nexMaxX) * 0.5d;
						extentsX = (nexMaxX - newMinX) * -0.5d;
					}
				}
			}
//		}
	}

	public void Encapsulate (Boundsd bounds)
	{
		Encapsulate(bounds.min);
		Encapsulate(bounds.max);
	}

	public void Expand (double amount)
	{
		//	Absolute value or minimum zero? Tough to say.
		m_extents.x = m_extents.x + 0.5 * amount;
		m_extents.y = m_extents.y + 0.5 * amount;
		m_extents.z = m_extents.z + 0.5 * amount;
	}

	public void Expand (Vector3d amounts)
	{
		m_extents.x = m_extents.x + 0.5 * amounts.x;
		m_extents.y = m_extents.y + 0.5 * amounts.y;
		m_extents.z = m_extents.z + 0.5 * amounts.z;
	}

	public bool IntersectRay (Rayd ray)
	{
//		if(Contains (ray.origin))
//			return true;
		double garbage;
		return IntersectRay(ray, out garbage);
	}

	public bool IntersectRay (Rayd ray, out double distance)
	{
		//	I could do this by creating six planes, and testing against all of them, but I only need three planes at most
		//	But I can do even better than that. Here:

		distance = 0d;
		Vector3d minHit, maxHit;
		if(!IntersectRayOnAxis(out minHit.x, out maxHit.x, ray.origin.x, ray.direction.x, center.x, extents.x))
			return false;
		if(!IntersectRayOnAxis(out minHit.y, out maxHit.y, ray.origin.y, ray.direction.y, center.y, extents.y))
			return false;
		if(!IntersectRayOnAxis(out minHit.z, out maxHit.z, ray.origin.z, ray.direction.z, center.z, extents.z))
			return false;
		double minDistance = Mathd.Max (minHit.x, minHit.y, minHit.z);
		double maxDistance = Mathd.Min (maxHit.x, maxHit.y, maxHit.z);
		if(minDistance > maxDistance) {
			return false;
		}
//		Debug.Log ("Distance: " + minDistance + " - " + maxDistance + "\n" + this + "  " + ray);
		if(minDistance <= 0d && maxDistance >= 0d){
//			distance = 0d;
			distance = minDistance;
			return true;
		} else if(minDistance > 0d){
			distance = minDistance;
			return true;
		} else {
			distance = maxDistance;
			distance = 0d;	//	Match original behavior for testing
			return false;
		}
	}

	bool IntersectRayOnAxis (out double minDistanceToHit, out double maxDistanceToHit, 
	                         double origin, double slope, 
	                         double center, double extents)
	{
		double originDiff = center - origin;
		if(slope == 0d){
			if(Mathd.Abs(originDiff) <= extents){
				minDistanceToHit = double.NegativeInfinity;
				maxDistanceToHit = double.PositiveInfinity;
				return true;
			} else {
				minDistanceToHit = double.PositiveInfinity;
				maxDistanceToHit = double.NegativeInfinity;
				return false;
			}
		} else {
			double inverseSlope = 1d / slope;
			double hit0 = (originDiff + extents) * inverseSlope;
			double hit1 = (originDiff - extents) * inverseSlope;
			if(hit0 < hit1){
				minDistanceToHit = hit0;
				maxDistanceToHit = hit1;
			} else {
				minDistanceToHit = hit1;
				maxDistanceToHit = hit0;
			}
			return true;
		}
	}

	public bool Intersects (Boundsd bounds)
	{
		//	Oh wait. This is easy.
		Vector3d diff = (bounds.center - this.center);
		Vector3d totalExtents = bounds.extents + this.extents;
		return Mathd.Abs(diff.x) < totalExtents.x
				&& Mathd.Abs(diff.y) < totalExtents.y
				&& Mathd.Abs(diff.z) < totalExtents.z;
	}

	public void SetMinMax (Vector3d min, Vector3d max)
	{
		this.center = (min + max) * 0.5;
		this.m_extents = (max - min) * 0.5;
	}

	public double SqrDistance (Vector3d point)
	{
		Vector3d diff = point - this.center;
		diff.x = Math.Max(0d, Math.Abs(diff.x) - m_extents.x);
		diff.y = Math.Max(0d, Math.Abs(diff.y) - m_extents.y);
		diff.z = Math.Max(0d, Math.Abs(diff.z) - m_extents.z);

		return diff.sqrMagnitude;
	}

	public double SqrDistance (Vector3d point, out Vector3d diff)
	{
		diff = point - this.center;
		diff.x = Math.Max(0d, Math.Abs(diff.x) - m_extents.x);
		diff.y = Math.Max(0d, Math.Abs(diff.y) - m_extents.y);
		diff.z = Math.Max(0d, Math.Abs(diff.z) - m_extents.z);
		
		return diff.sqrMagnitude;
	}

	public double SqrDistanceNotNaively (Vector3d point)
	{
		Vector3d diff = point - this.center;
		//	Subtract the extents but cap at zero.
		//	Do it naively. Assume extents are always positive.
		//	Crap. No, I kind of need absolute values all the way around.
		diff.x = Math.Max(0d, Math.Abs(diff.x) - Math.Abs (m_extents.x));
		diff.y = Math.Max(0d, Math.Abs(diff.y) - Math.Abs (m_extents.y));
		diff.z = Math.Max(0d, Math.Abs(diff.z) - Math.Abs (m_extents.z));
		
		return diff.sqrMagnitude;
	}

	public double SqrDistanceNotNaively (Vector3d point, out Vector3d diff)
	{
		diff = point - this.center;
		//	Subtract the extents but cap at zero.
		//	Do it naively. Assume extents are always positive.
		//	Crap. No, I kind of need absolute values all the way around.
		diff.x = Math.Max(0d, Math.Abs(diff.x) - Math.Abs (m_extents.x));
		diff.y = Math.Max(0d, Math.Abs(diff.y) - Math.Abs (m_extents.y));
		diff.z = Math.Max(0d, Math.Abs(diff.z) - Math.Abs (m_extents.z));
		
		return diff.sqrMagnitude;
	}
}

