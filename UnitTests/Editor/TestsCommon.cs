using System;
using UnityEngine;
using System.Collections.Generic;
using NUnit.Framework;
using System.Text;



namespace UnityTest
{
	public class TestsCommon
	{
		public const float floatPrecisionF = 1.19209289550781E-07f;
		public const double floatPrecision = 1.19209289550781E-07;
		public const double doublePrecision = 2.22044604925031E-16;
		public const float floatBumpUp = 1.00000012f;
		public const float floatBumpDown= 0.99999994f;	//	It's worth noting that the bump down is an extra factor of 2 smaller
		public const double doubleBumpUp = 1.0000000000000002d;
		public const double doubleBumpDown = 0.99999999999999989d;
		
		public const string fullNumberFormat = "G7";
		public const string fullBoundsFormatOneLine = "{center}\u00B1{extents}|{min}to{max}";
		public const string fullBoundsFormat = "{center}\u00B1{extents}\n{min}to{max}";


		public TestsCommon ()
		{
			TestConstEqual("floatBumpUp", floatBumpUp, 1f + floatPrecisionF);
			TestConstEqual("floatBumpDown", floatBumpDown, 1f - 0.5f * floatPrecisionF);
			TestConstEqual("doubleBumpUp", doubleBumpUp, 1d + doublePrecision);
			TestConstEqual("doubleBumpDown", doubleBumpDown, 1d - 0.5d * doublePrecision);
		}

		void TestConstEqual (string name, double actual, double expected)
		{
			if(actual != expected){
				Debug.LogError ("Constant value " + name + " does not match expected value. " +
					"Actual=" + actual.ToString("G30") + "  Expected=" + expected.ToString("G30"));
			}
		}

		protected virtual void OnAssertSimilarFailure(double f, double d, string valueName, double toleranceBasedOn, double difference, double tolerance){
			string debugString = "f=" + f + ", \td=" + d + "\n\t"
				+ "diff=" + difference + ", \ttolerance=" + tolerance + "\n\t"
					+ "factors= " + (tolerance == 0.0d ? "&#x221E" : (difference / tolerance).ToString())
					+ " \t/\t " + (difference == 0.0d ? "&#x221E" : (tolerance / difference).ToString());
			Assert.Fail(valueName + " outside of tolerance\n\t" + debugString, f, d);
		}
		
		protected virtual void AssertSimilar(float f, double d){
			AssertSimilar(f, d, "Values");
		}
		
		protected virtual void AssertSimilar(float f, double d, double toleranceBasedOn){
			AssertSimilar(f, d, "Values", toleranceBasedOn);
		}
		protected virtual void AssertSimilar(float f, double d, String valueName){
			AssertSimilar(f, d, valueName, 1d);
		}
		protected virtual void AssertSimilar(double f, double d, string valueName, double toleranceBasedOn){
			double difference = Mathd.Abs(d - f);
			double tolerance = Mathd.Abs(toleranceBasedOn * floatPrecision);
			if(difference > tolerance){
				OnAssertSimilarFailure(f, d, valueName, toleranceBasedOn, difference, tolerance);
			} else {
				//				if(difference != 0){
				//					double passFactor = Mathd.Min ((tolerance / floatPrecision), (tolerance / difference));
				//					if(passFactor >= 10.0d){
				//						Debug.Log (valueName + " inside of tolerance (" + debugString + ")");
				//					}
				//				}
				//				Assert.Pass(valueName + " inside of tolerance", f, d);
			}
		}

		//	Vectors
		protected virtual void AssertSimilar(Vector3 f, Vector3 d){
			AssertSimilar(f, d, 1);
		}
		protected virtual void AssertSimilar(Vector3 f, Vector3d d){
			AssertSimilar(f, d, 1);
		}
		protected virtual void AssertSimilar(Vector3d f, Vector3d d){
			AssertSimilar(f, d, 1);
		}
		
		protected virtual void AssertSimilar(Vector3 f, Vector3 d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		protected virtual void AssertSimilar(Vector3 f, Vector3d d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		protected virtual void AssertSimilar(Vector3d f, Vector3d d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		
		protected virtual void AssertSimilar(Vector3 f, Vector3 d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		protected virtual void AssertSimilar(Vector3 f, Vector3d d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		protected virtual void AssertSimilar(Vector3d f, Vector3d d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		
		protected virtual void AssertSimilar(Vector3 f, Vector3 d, double toleranceBasedOn, string additionalInfo){
			AssertSimilar(new Vector3d(f), new Vector3d(d), toleranceBasedOn, additionalInfo);
		}
		
		protected virtual void AssertSimilar(Vector3 f, Vector3d d, double toleranceBasedOn, string additionalInfo){
			AssertSimilar(new Vector3d(f), d, toleranceBasedOn, additionalInfo);
		}
		
		
		protected virtual void AssertSimilar(Vector3d f, Vector3d d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
		}

		//	Quaternions
		protected virtual void AssertSimilar(Quaternion f, Quaternion d){
			AssertSimilar(f, d, 1);
		}
		protected virtual void AssertSimilar(Quaternion f, Quaterniond d){
			AssertSimilar(f, d, 1);
		}
		protected virtual void AssertSimilar(Quaterniond f, Quaterniond d){
			AssertSimilar(f, d, 1);
		}

		protected virtual void AssertSimilar(Quaternion f, Quaternion d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		protected virtual void AssertSimilar(Quaternion f, Quaterniond d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		protected virtual void AssertSimilar(Quaterniond f, Quaterniond d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}

		protected virtual void AssertSimilar(Quaternion f, Quaternion d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		protected virtual void AssertSimilar(Quaternion f, Quaterniond d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		protected virtual void AssertSimilar(Quaterniond f, Quaterniond d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		
		protected virtual void AssertSimilar(Quaternion f, Quaternion d, double toleranceBasedOn, string additionalInfo){
			AssertSimilar(new Quaterniond(f), new Quaterniond(d), toleranceBasedOn, additionalInfo);
		}
		protected virtual void AssertSimilar(Quaternion f, Quaterniond d, double toleranceBasedOn, string additionalInfo){
			AssertSimilar(new Quaterniond(f), d, toleranceBasedOn, additionalInfo);
		}
		protected virtual void AssertSimilar(Quaterniond f, Quaterniond d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "  d: " + d.ToString("G5") + "\n";
			Quaterniond dneg = Negative(d);
			if(GetSumOfSquares(Subtract(f, d)) > GetSumOfSquares(Subtract(f, dneg))){
				d = dneg;
			}
			AssertSimilar (f.x, d.x, inputs + "x", toleranceBasedOn);
			AssertSimilar (f.y, d.y, inputs + "y", toleranceBasedOn);
			AssertSimilar (f.z, d.z, inputs + "z", toleranceBasedOn);
			AssertSimilar (f.w, d.w, inputs + "w", toleranceBasedOn);
		}

		//	Bounds
		protected virtual void AssertSimilar(Bounds f, Bounds d){
			AssertSimilar(f, d, 1);
		}
		protected virtual void AssertSimilar(Bounds f, Boundsd d){
			AssertSimilar(f, d, 1);
		}
		protected virtual void AssertSimilar(Boundsd f, Boundsd d){
			AssertSimilar(f, d, 1);
		}
		
		protected virtual void AssertSimilar(Bounds f, Bounds d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		protected virtual void AssertSimilar(Bounds f, Boundsd d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		protected virtual void AssertSimilar(Boundsd f, Boundsd d, double toleranceBasedOn){
			AssertSimilar(f, d, toleranceBasedOn, "");
		}
		
		protected virtual void AssertSimilar(Bounds f, Bounds d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		protected virtual void AssertSimilar(Bounds f, Boundsd d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		protected virtual void AssertSimilar(Boundsd f, Boundsd d, string additionalInfo){
			AssertSimilar(f, d, 1d, additionalInfo);
		}
		
		protected virtual void AssertSimilar(Bounds f, Bounds d, double toleranceBasedOn, string additionalInfo){
			AssertSimilar(new Boundsd(f), new Boundsd(d), toleranceBasedOn, additionalInfo);
		}
		protected virtual void AssertSimilar(Bounds f, Boundsd d, double toleranceBasedOn, string additionalInfo){
			AssertSimilar(new Boundsd(f), d, toleranceBasedOn, additionalInfo);
		}
		protected virtual void AssertSimilar(Boundsd f, Boundsd d, double toleranceBasedOn, string additionalInfo){
			string inputs = (additionalInfo.Length > 1 ? additionalInfo + "\n" : "") + "f: " + f.ToString("G5") + "\nd: " + d.ToString("G5") + "\ntesting ";
			AssertSimilar (f.center, d.center, toleranceBasedOn, inputs + "center");
			AssertSimilar (f.extents, d.extents, toleranceBasedOn, inputs + "extents");
			AssertSimilar (f.min, d.min, toleranceBasedOn, inputs + "min");
			AssertSimilar (f.max, d.max, toleranceBasedOn, inputs + "max");
			AssertSimilar (f.size, d.size, toleranceBasedOn);
		}

		public static Quaterniond Negative(Quaterniond q){
			return new Quaterniond(-q.x, -q.y, -q.z, -q.w);
		}
		public static Quaterniond Subtract(Quaternion a, Quaterniond b){
			return new Quaterniond(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
		}
		public static Quaterniond Subtract(Quaterniond a, Quaterniond b){
			return new Quaterniond(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
		}

		public static void Normalize(ref Quaternion q){
			float sum = 0;
			for(int i=0; i<4; ++i){
				sum += q[i] * q[i];
			}
			float inverseMagnitude = 1f / Mathf.Sqrt(sum);
			for(int i=0; i<4; ++i){
				q[i] *= inverseMagnitude;
			}
		}
		
		public static bool IsNormalized(Quaternion q, double tolerance){
			return Math.Abs(1d - GetSumOfSquares(q)) < tolerance;
		}
		public static bool IsNormalized(Quaterniond q, double tolerance){
			return Math.Abs(1d - GetSumOfSquares(q)) < tolerance;
		}
		public static double GetSumOfSquares(Quaternion q){
			double sum = 0;
			for(int i=0; i<4; ++i){
				sum += ((double)q[i]) * q[i];
			}
			return sum;
		}
		public static double GetSumOfSquares(Quaterniond q){
			double sum = 0;
			for(int i=0; i<4; ++i){
				sum += q[i] * q[i];
			}
			return sum;
		}

		public static Quaternion GenerateRandomQuaternion(System.Random rand){
			Quaternion q = new Quaternion();
			for(int i=0; i<4; ++i){
				q[i] = (float)(2.0 * rand.NextDouble() - 1.0);	//	I want a range of negative and positive values.
			}
			Normalize (ref q);	//	Not sure if this will actually be valid or not. If it produces errors, it'll be funny at least.
			return q;
		}
		public static Vector3 GenerateRandomVector3(System.Random rand){
			Vector3 v3 = new Vector3();
			for(int i=0; i<3; ++i){
				//	I want a range of negative and positive values, with a hilariously large range of magnitudes.
				v3[i] = (float)(1000.0 * (0.1 + rand.NextDouble()) * (2.0 * rand.NextDouble() - 1.0));	
			}
			return v3;
		}

		
		public static Bounds GenerateRandomBounds (System.Random rand)
		{
			return new Bounds(GenerateRandomVector3(rand), GenerateRandomVector3(rand));
		}


		public static Vector3[] RollPointsWithinBounds (System.Random rand, Bounds bounds, int count, float extentRatio)
		{
			bounds.extents = bounds.extents * extentRatio;
			return RollPointsWithinBounds(rand, bounds, count);
		}
		public static Vector3[] RollPointsWithinBounds (System.Random rand, Bounds bounds, int count)
		{
			Vector3[] points = new Vector3[count];
			for(int i=0; i<count; ++i){
				points[i] = RollPointWithinBounds(rand, bounds);
			}
			return points;
		}
		
		public static Vector3 RollPointWithinBounds (System.Random rand, Bounds bounds)
		{
			return new Vector3(bounds.center.x + RollFloatBetween(rand, -bounds.extents.x, bounds.extents.x),
			                   bounds.center.y + RollFloatBetween(rand, -bounds.extents.y, bounds.extents.y),
			                   bounds.center.z + RollFloatBetween(rand, -bounds.extents.z, bounds.extents.z));
		}
		
		public static float RollFloatBetween (System.Random rand, double min, double max)
		{
			return (float)RollDoubleBetween(rand, min, max);
		}
		public static double RollDoubleBetween (System.Random rand, double min, double max)
		{
			return min + rand.NextDouble() * (max - min);
		}

		public static Vector3[] RollPointsOutsideBounds (System.Random rand, Bounds bounds, int count, float extentRatio)
		{
			bounds.extents = bounds.extents * extentRatio;
			return RollPointsOutsideBounds(rand, bounds, count);
		}

		public static Vector3[] RollPointsOutsideBounds (System.Random rand, Bounds bounds, int count)
		{
			Vector3[] points = new Vector3[count];
			for(int i=0; i<count; ++i){
				points[i] = RollPointOutsideBounds(rand, bounds);
			}
			return points;
		}
		
		public static Vector3 RollPointOutsideBounds (System.Random rand, Bounds bounds)
		{
			//	Now, some of the coordinates can be inside the bounds. They just can't ALL be inside.
			//	I want to generate all possibilities, but with some skewing.
			int inPatternNum = rand.Next() % 7;	//	This excludes the binary value 111, while allowing every other combination.
			bool[] pattern = {inPatternNum % 2 == 1, (inPatternNum / 2) % 2 == 1,(inPatternNum / 4) % 2 == 1};
			Vector3 point = new Vector3();
			for(int i=0; i<3; ++i){
				if(pattern[i]){
					point[i] = bounds.center[i] + RollFloatBetween(rand, -bounds.extents[i], bounds.extents[i]);
				} else {
					point[i] = RollOutsideBoundsAxis(rand, bounds.center[i], Mathf.Abs(bounds.extents[i]));
				}
			}
			return point;
		}

		public static float RollOutsideBoundsAxis (System.Random rand, float center, float extent)
		{
			float roll = (float)(rand.NextDouble() * extent);

			if(rand.Next() % 2 == 0){
				//	Below minimum
				float value = center - extent - roll;
				if(center - value <= extent){
					float firstValue = value;
					value = bumpFloor(center - extent);
					if(center - value <= extent){
						Debug.LogError("Problem rolling outside bounds axis (center=" + center 
						           + ";extent=" + extent +";roll=" + roll 
						           + ";firstValue=" + firstValue + ";value=" + value + ")");
					} else {
						Debug.LogWarning("Difficulty rolling outside bounds axis (center=" + center 
						               + ";extent=" + extent +";roll=" + roll 
						               + ";firstValue=" + firstValue + ";value=" + value + ")");
					}
				}
				return value;
			} else {
				//	Above maximum
				float value = center + extent + roll;
				if(value - center <= extent){
					float firstValue = value;
					value = bumpCeiling(center + extent);
					if(value - center <= extent){
						Debug.LogError("Problem rolling outside bounds axis (center=" + center 
						               + ";extent=" + extent +";roll=" + roll 
						               + ";firstValue=" + firstValue + ";value=" + value + ")");
					} else {
						Debug.LogWarning("Difficulty rolling outside bounds axis (center=" + center 
						                 + ";extent=" + extent +";roll=" + roll 
						                 + ";firstValue=" + firstValue + ";value=" + value + ")");
					}
				}
				return value;
			}
		}

		public static float bumpFloor(float value){
			if(value == 0f){
				return -float.MinValue;
			} else if(value < 0f){
				return bumpUp(value);
			} else {
				return bumpDown(value);
			}
		}
		public static float bumpCeiling(float value){
			if(value == 0f){
				return float.MinValue;
			} else if(value < 0f){
				return bumpDown(value);
			} else {
				return bumpUp(value);
			}
		}
		public static float bumpUp(float value){
			if(MathUtil.IsNegativeZero(value)){
				return -float.MinValue;
			} else if(value == 0f) {
				return float.MinValue;
			} else if(float.IsNaN(value)){
				return value;
			} else if(float.IsInfinity(value)){
				return value;
			} else if(value == float.MaxValue) {
				return float.PositiveInfinity;
			} else if(value == -float.MaxValue) {
				return float.NegativeInfinity;
			} else {
				float result = value * floatBumpUp;
				if(result == value){	//	Uh oh! The single bump didn't do it!
					float bumpDelta = floatPrecisionF;
					int multiple = 1;
					while(result == value){
						bumpDelta *= 2;
						multiple *= 2;
						result = value * (1f + bumpDelta);
						if(bumpDelta > value){
							throw new Exception("Bump Delta reached " + bumpDelta + " at 2^" + multiple + 
							               " with result=" + result + " and value=" + value + " before emergency termination");
							//	WTF. In fact, a while loop shouldn't be necessary. Just multiplying by two will work just fine.
//							return result;
						}
					}
					Debug.Log ("Single bump didn't work on " + value.ToString("G30") + ", required 2^" 
					           + multiple + " bump factors to resolve to: " + result.ToString("G30") + " (bumpDelta=" + bumpDelta + ")");
				}
				return result;
			}
		}
		public static float bumpDown(float value){
			if(value == 0f){
				return 0f;
			} else if(float.IsNaN(value)){
				return value;
			} else if(float.IsInfinity(value)){
				return value;
			} else if(value == float.MinValue) {
				return 0f;
			} else if(value == -float.MinValue) {
				return -0.0f;
			} else {
				float result = value * floatBumpDown;
				if(result == value){	//	Uh oh! The single bump didn't do it!
					float bumpDelta = -floatPrecisionF;
					float multiple = 1;
					while(result == value){
						bumpDelta *= 2;
						multiple++;
						result = value * (1f + bumpDelta);
						if(multiple >= 27){
							//	WTF. In fact, a while loop shouldn't be necessary. Just multiplying by two should always be enough.
							throw new Exception("Bump Delta reached " + bumpDelta + " at 2^" + multiple + 
							               " with result=" + result + " and value=" + value + " before emergency termination");
//							return result;
						}
					}
					Debug.Log ("Single bump didn't work on " + value.ToString("G30") + ", required 2^" 
					           + multiple + " bump factors to resolve to: " + result.ToString("G30") + " (bumpDelta=" + bumpDelta + ")");
				}
				return result;
			}
		}

		//	ToString shortcuts for getting maximum data
		
		public static string Str(Vector3 v3){
			return v3.ToString(fullNumberFormat);
		}
		public static string Str(Vector3d v3){
			return v3.ToString(fullNumberFormat);
		}
		public static string Str(Quaternion q){
			return q.ToString(fullNumberFormat);
		}
		public static string Str(Quaterniond q){
			return q.ToString(fullNumberFormat);
		}
		public static string Str(Bounds b){
			return Boundsd.ToString(b, fullBoundsFormatOneLine, fullNumberFormat);
		}
		public static string Str(Boundsd b){
			return b.ToString(fullBoundsFormatOneLine, fullNumberFormat);
		}
		public static string Str(Ray ray){
			return ray.ToString(fullNumberFormat);
		}
		public static string Str(Rayd ray){
			return ray.ToString(fullNumberFormat);
		}
		public static string StrMultiline(Bounds b){
			return StrMultiline(b, 1);
		}
		public static string StrMultiline(Boundsd b){
			return StrMultiline(b, 1);
		}
		public static string StrMultiline(Bounds b, int tabIndentsPerLine){
			string newlineReplacement = "\n" + StringUtil.Repeat("\t", tabIndentsPerLine);
			return Boundsd.ToString(b, fullBoundsFormat.Replace("\n", newlineReplacement), fullNumberFormat);
		}
		public static string StrMultiline(Boundsd b, int tabIndentsPerLine){
			string newlineReplacement = "\n" + StringUtil.Repeat("\t", tabIndentsPerLine);
			return b.ToString(fullBoundsFormat.Replace("\n", newlineReplacement), fullNumberFormat);
		}
	}
}

