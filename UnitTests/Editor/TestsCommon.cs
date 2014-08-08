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

		public NumberRollProfile worldCoordinateRoller, valueRoller, quaternionRoller, 
				eulerRoller, boundsCenterRoller, boundsExtentsRoller,
				rayOriginRoller, rayDirectionRoller;

		public TestsCommon ()
		{
			TestConstants();
			InitializeRollProfiles();
		}

		void InitializeRollProfiles ()
		{
			valueRoller = new NumberRollProfile();
			worldCoordinateRoller = new NumberRollProfile(Rolls.mult(Rolls.r(10,1100),Rolls.r(-1,1)));
			quaternionRoller = new NumberRollProfile(-1f, 1f);
			eulerRoller = new NumberRollProfile(Rolls.r (0f, 360f))
				.AddOverride(0, new NumberRollProfile(Rolls.a(Rolls.r(0,90),Rolls.s(0, 270))));
			boundsCenterRoller = worldCoordinateRoller.Clone();
			boundsExtentsRoller = new NumberRollProfile(Rolls.mult(Rolls.r(10,1100),Rolls.r(-1,1)));
			rayOriginRoller = worldCoordinateRoller.Clone();
			rayDirectionRoller = new NumberRollProfile(-1f, 1f);
		}

		//	Object generation for test sets
		public Quaternion GenerateRandomQuaternion(System.Random rand){
			return GenerateRandomQuaternion(rand, quaternionRoller);
		}
		public Vector3 GenerateRandomVector3(System.Random rand){
			return GenerateRandomVector3(rand, worldCoordinateRoller);
		}
		public Vector3 GenerateRandomEulerCoordinates(System.Random rand){
			return GenerateRandomVector3(rand, eulerRoller);
		}
		public Bounds GenerateRandomBounds(System.Random rand){
			return GenerateRandomBounds(rand, boundsCenterRoller, boundsExtentsRoller);
		}
		public Ray GenerateRandomRay(System.Random rand){
			return GenerateRandomRay(rand, rayOriginRoller, rayDirectionRoller);
		}

		public static Quaternion GenerateRandomQuaternion(System.Random rand, NumberRollProfile roller){
			Quaternion q = new Quaternion();
			for(int i=0; i<4; ++i){
				q[i] = roller.Roll(i, rand);
			}
			Normalize (ref q);
			return q;
		}
		public static Vector3 GenerateRandomVector3(System.Random rand, NumberRollProfile roller){
			Vector3 v3 = new Vector3();
			for(int i=0; i<3; ++i){
				//	I want a range of negative and positive values, with a hilariously large range of magnitudes.
				v3[i] = roller.Roll(i, rand);
			}
			return v3;
		}
		
		public static Bounds GenerateRandomBounds (System.Random rand, NumberRollProfile centerRoller, NumberRollProfile extentsRoller)
		{
			return new Bounds(GenerateRandomVector3(rand, centerRoller), GenerateRandomVector3(rand, extentsRoller));
		}
		public static Ray GenerateRandomRay(System.Random rand, NumberRollProfile originRoller, NumberRollProfile directionRoller)
		{
			return new Ray(GenerateRandomVector3(rand, originRoller), MiscUtil.RollRandomPointOnSphereF(rand));
		}

		//	Old versions. Don't want to lose them until I'm sure I don't need them.
		public static Quaternion GenerateRandomQuaternionOld(System.Random rand){
			Quaternion q = new Quaternion();
			for(int i=0; i<4; ++i){
				q[i] = (float)(2.0 * rand.NextDouble() - 1.0);	//	I want a range of negative and positive values.
			}
			Normalize (ref q);	//	Not sure if this will actually be valid or not. If it produces errors, it'll be funny at least.
			return q;
		}
		public static Vector3 GenerateRandomVector3Old(System.Random rand){
			Vector3 v3 = new Vector3();
			for(int i=0; i<3; ++i){
				//	I want a range of negative and positive values, with a hilariously large range of magnitudes.
				v3[i] = (float)(1000.0 * (0.1 + rand.NextDouble()) * (2.0 * rand.NextDouble() - 1.0));	
			}
			return v3;
		}
		
		public static Bounds GenerateRandomBoundsOld (System.Random rand)
		{
			return new Bounds(GenerateRandomVector3Old(rand), GenerateRandomVector3Old(rand));
		}



		void TestConstants ()
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
			if(double.IsNaN(f) || double.IsNaN(f) || double.IsInfinity(f) || double.IsInfinity(d)){
				if(!f.Equals(d)){
					OnAssertSimilarFailure(f, d, valueName, toleranceBasedOn, 0d, 0d);
				}
			} else {

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
				if(!(float.IsNaN(q[i]) || float.IsInfinity(q[i]))){
					sum += q[i] * q[i];
				}
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
			if(float.IsInfinity(extent)){
				//	Not possible to pick a location outside the bounds axis
				//	Unless I go with NaN. I'm not sure that makes sense.
				return center;
			}
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

		public static bool ContainsNaN(Bounds b){
			return ContainsNaN(b.center) || ContainsNaN(b.extents);
		}
		public static bool ContainsNaN(Vector3 v){
			return float.IsNaN(v.x) || float.IsNaN(v.y) || float.IsNaN(v.z);
		}
		public static bool ContainsNaN(Quaternion q){
			return float.IsNaN(q.x) || float.IsNaN(q.y) || float.IsNaN(q.z) || float.IsNaN(q.w);
		}
		public static bool ContainsNaN(Ray r){
			return ContainsNaN(r.origin) || ContainsNaN(r.direction);
		}
		public static bool ContainsNaN(Vector2 v){
			return float.IsNaN(v.x) || float.IsNaN(v.y);
		}
		public static bool ContainsNaN(float f){
			return float.IsNaN(f);
		}
		public static bool ContainsNaN(params float[] values){
			foreach(float value in values){
				if(float.IsNaN(value))
					return true;
			}
			return false;
		}

		public static bool ContainsNaN(Boundsd b){
			return ContainsNaN(b.center) || ContainsNaN(b.extents);
		}
		public static bool ContainsNaN(Vector3d v){
			return double.IsNaN(v.x) || double.IsNaN(v.y) || double.IsNaN(v.z);
		}
		public static bool ContainsNaN(Quaterniond q){
			return double.IsNaN(q.x) || double.IsNaN(q.y) || double.IsNaN(q.z) || double.IsNaN(q.w);
		}
		public static bool ContainsNaN(Rayd r){
			return ContainsNaN(r.origin) || ContainsNaN(r.direction);
		}
		public static bool ContainsNaN(Vector2d v){
			return double.IsNaN(v.x) || double.IsNaN(v.y);
		}
		public static bool ContainsNaN(double d){
			return double.IsNaN(d);
		}
		public static bool ContainsNaN(params double[] values){
			foreach(double value in values){
				if(double.IsNaN(value))
					return true;
			}
			return false;
		}

		public static bool ContainsNaN (object obj)
		{
			if(obj == null){
				return false;
			}
			if(obj is Bounds){
				return ContainsNaN((Bounds)obj);
			} else if(obj is Vector3){
				return ContainsNaN((Vector3)obj);
			} else if(obj is Quaternion){
				return ContainsNaN((Quaternion)obj);
			} else if(obj is Ray){
				return ContainsNaN((Ray)obj);
			} else if(obj is Vector2){
				return ContainsNaN((Vector2)obj);
			} else if(obj is Single){
				return ContainsNaN((float)obj);
			}

			if(obj is Boundsd){
				return ContainsNaN((Boundsd)obj);
			} else if(obj is Vector3d){
				return ContainsNaN((Vector3d)obj);
			} else if(obj is Quaterniond){
				return ContainsNaN((Quaterniond)obj);
			} else if(obj is Rayd){
				return ContainsNaN((Rayd)obj);
			} else if(obj is Vector2d){
				return ContainsNaN((Vector2d)obj);
			} else if(obj is Double){
				return ContainsNaN((double)obj);
			}

			Debug.LogError("ContainsNaN couldn't resolve: Unknown object type: " + obj.GetType() + " (obj='" + obj + "')");
			return false;
		}
		public static bool ContainsNaN (params object[] objects)
		{
			foreach(object obj in objects){
				if(ContainsNaN(obj))
					return true;
			}
			return false;
		}

		public static bool ContainsPositiveInfinity(Bounds b){
			return ContainsPositiveInfinity(b.center) || ContainsPositiveInfinity(b.extents);
		}
		public static bool ContainsPositiveInfinity(Vector3 v){
			return float.IsPositiveInfinity(v.x) || float.IsPositiveInfinity(v.y) || float.IsPositiveInfinity(v.z);
		}
		public static bool ContainsPositiveInfinity(Quaternion q){
			return float.IsPositiveInfinity(q.x) || float.IsPositiveInfinity(q.y) || float.IsPositiveInfinity(q.z) || float.IsPositiveInfinity(q.w);
		}
		public static bool ContainsPositiveInfinity(Ray r){
			return ContainsPositiveInfinity(r.origin) || ContainsPositiveInfinity(r.direction);
		}
		public static bool ContainsPositiveInfinity(Vector2 v){
			return float.IsPositiveInfinity(v.x) || float.IsPositiveInfinity(v.y);
		}
		public static bool ContainsPositiveInfinity(float f){
			return float.IsPositiveInfinity(f);
		}
		public static bool ContainsPositiveInfinity(params float[] values){
			foreach(float value in values){
				if(float.IsPositiveInfinity(value))
					return true;
			}
			return false;
		}
		
		public static bool ContainsPositiveInfinity(Boundsd b){
			return ContainsPositiveInfinity(b.center) || ContainsPositiveInfinity(b.extents);
		}
		public static bool ContainsPositiveInfinity(Vector3d v){
			return double.IsPositiveInfinity(v.x) || double.IsPositiveInfinity(v.y) || double.IsPositiveInfinity(v.z);
		}
		public static bool ContainsPositiveInfinity(Quaterniond q){
			return double.IsPositiveInfinity(q.x) || double.IsPositiveInfinity(q.y) || double.IsPositiveInfinity(q.z) || double.IsPositiveInfinity(q.w);
		}
		public static bool ContainsPositiveInfinity(Rayd r){
			return ContainsPositiveInfinity(r.origin) || ContainsPositiveInfinity(r.direction);
		}
		public static bool ContainsPositiveInfinity(Vector2d v){
			return double.IsPositiveInfinity(v.x) || double.IsPositiveInfinity(v.y);
		}
		public static bool ContainsPositiveInfinity(double d){
			return double.IsPositiveInfinity(d);
		}
		public static bool ContainsPositiveInfinity(params double[] values){
			foreach(double value in values){
				if(double.IsPositiveInfinity(value))
					return true;
			}
			return false;
		}
		
		public static bool ContainsPositiveInfinity (object obj)
		{
			if(obj == null){
				return false;
			}
			if(obj is Bounds){
				return ContainsPositiveInfinity((Bounds)obj);
			} else if(obj is Vector3){
				return ContainsPositiveInfinity((Vector3)obj);
			} else if(obj is Quaternion){
				return ContainsPositiveInfinity((Quaternion)obj);
			} else if(obj is Ray){
				return ContainsPositiveInfinity((Ray)obj);
			} else if(obj is Vector2){
				return ContainsPositiveInfinity((Vector2)obj);
			} else if(obj is Single){
				return ContainsPositiveInfinity((float)obj);
			}
			
			if(obj is Boundsd){
				return ContainsPositiveInfinity((Boundsd)obj);
			} else if(obj is Vector3d){
				return ContainsPositiveInfinity((Vector3d)obj);
			} else if(obj is Quaterniond){
				return ContainsPositiveInfinity((Quaterniond)obj);
			} else if(obj is Rayd){
				return ContainsPositiveInfinity((Rayd)obj);
			} else if(obj is Vector2d){
				return ContainsPositiveInfinity((Vector2d)obj);
			} else if(obj is Double){
				return ContainsPositiveInfinity((double)obj);
			}
			
			Debug.LogError("ContainsPositiveInfinity couldn't resolve: Unknown object type: " + obj.GetType() + " (obj='" + obj + "')");
			return false;
		}
		public static bool ContainsPositiveInfinity (params object[] objects)
		{
			foreach(object obj in objects){
				if(ContainsPositiveInfinity(obj))
					return true;
			}
			return false;
		}

		public static bool ContainsNegativeInfinity(Bounds b){
			return ContainsNegativeInfinity(b.center) || ContainsNegativeInfinity(b.extents);
		}
		public static bool ContainsNegativeInfinity(Vector3 v){
			return float.IsNegativeInfinity(v.x) || float.IsNegativeInfinity(v.y) || float.IsNegativeInfinity(v.z);
		}
		public static bool ContainsNegativeInfinity(Quaternion q){
			return float.IsNegativeInfinity(q.x) || float.IsNegativeInfinity(q.y) || float.IsNegativeInfinity(q.z) || float.IsNegativeInfinity(q.w);
		}
		public static bool ContainsNegativeInfinity(Ray r){
			return ContainsNegativeInfinity(r.origin) || ContainsNegativeInfinity(r.direction);
		}
		public static bool ContainsNegativeInfinity(Vector2 v){
			return float.IsNegativeInfinity(v.x) || float.IsNegativeInfinity(v.y);
		}
		public static bool ContainsNegativeInfinity(float f){
			return float.IsNegativeInfinity(f);
		}
		public static bool ContainsNegativeInfinity(params float[] values){
			foreach(float value in values){
				if(float.IsNegativeInfinity(value))
					return true;
			}
			return false;
		}
		
		public static bool ContainsNegativeInfinity(Boundsd b){
			return ContainsNegativeInfinity(b.center) || ContainsNegativeInfinity(b.extents);
		}
		public static bool ContainsNegativeInfinity(Vector3d v){
			return double.IsNegativeInfinity(v.x) || double.IsNegativeInfinity(v.y) || double.IsNegativeInfinity(v.z);
		}
		public static bool ContainsNegativeInfinity(Quaterniond q){
			return double.IsNegativeInfinity(q.x) || double.IsNegativeInfinity(q.y) || double.IsNegativeInfinity(q.z) || double.IsNegativeInfinity(q.w);
		}
		public static bool ContainsNegativeInfinity(Rayd r){
			return ContainsNegativeInfinity(r.origin) || ContainsNegativeInfinity(r.direction);
		}
		public static bool ContainsNegativeInfinity(Vector2d v){
			return double.IsNegativeInfinity(v.x) || double.IsNegativeInfinity(v.y);
		}
		public static bool ContainsNegativeInfinity(double d){
			return double.IsNegativeInfinity(d);
		}
		public static bool ContainsNegativeInfinity(params double[] values){
			foreach(double value in values){
				if(double.IsNegativeInfinity(value))
					return true;
			}
			return false;
		}
		
		public static bool ContainsNegativeInfinity (object obj)
		{
			if(obj == null){
				return false;
			}
			if(obj is Bounds){
				return ContainsNegativeInfinity((Bounds)obj);
			} else if(obj is Vector3){
				return ContainsNegativeInfinity((Vector3)obj);
			} else if(obj is Quaternion){
				return ContainsNegativeInfinity((Quaternion)obj);
			} else if(obj is Ray){
				return ContainsNegativeInfinity((Ray)obj);
			} else if(obj is Vector2){
				return ContainsNegativeInfinity((Vector2)obj);
			} else if(obj is Single){
				return ContainsNegativeInfinity((float)obj);
			}
			
			if(obj is Boundsd){
				return ContainsNegativeInfinity((Boundsd)obj);
			} else if(obj is Vector3d){
				return ContainsNegativeInfinity((Vector3d)obj);
			} else if(obj is Quaterniond){
				return ContainsNegativeInfinity((Quaterniond)obj);
			} else if(obj is Rayd){
				return ContainsNegativeInfinity((Rayd)obj);
			} else if(obj is Vector2d){
				return ContainsNegativeInfinity((Vector2d)obj);
			} else if(obj is Double){
				return ContainsNegativeInfinity((double)obj);
			}
			
			Debug.LogError("ContainsNegativeInfinity couldn't resolve: Unknown object type: " + obj.GetType() + " (obj='" + obj + "')");
			return false;
		}
		public static bool ContainsNegativeInfinity (params object[] objects)
		{
			foreach(object obj in objects){
				if(ContainsNegativeInfinity(obj))
					return true;
			}
			return false;
		}
	}
}

