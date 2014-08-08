using System;
using UnityEngine;
using System.Collections.Generic;
using NUnit.Framework;
using System.Text;

namespace UnityTest {
	[TestFixture]
	[Category ("Bounds")]
	public class BoundsdTests : TestsCommon
	{
		const int numberOfTestItems = 30;
		const int expandedNumberOfTestItems = 60;
		const bool prenormalizeExtents = false;
		const bool useExtremeValuesInBoundsCenter = false;
		const bool useExtremeValuesInBoundsExtents = false;
		const bool useExtremeValuesInPoints = false;
		const bool useExtremeValuesInRayOrigins = false;
		const bool useExtremeValuesInRayDirections = false;
		const bool useExtremeValuesInRawValues = false;
		
		TestItemSet[] testItemSets;
		
		public struct TestItemSet {
			public Quaternion fq0, fq1;
			public Quaterniond dq0, dq1;
			public Vector3 fv0, fv1;
			public Vector3d dv0, dv1;
			public Bounds fb0, fb1;
			public Boundsd db0, db1;
			public float f0;
			public double d0;
			public Ray fray;
			public Rayd dray;

			public int seed;

			public override string ToString(){
				return ToString("G5");
			}
			public string ToString(string format){
				return "q0=" + dq0.ToString(format) + " q1=" + dq1.ToString(format) 
					+ " v0=" + dv0.ToString(format) + " v1=" + dv1.ToString(format) 
						+ " d0=" + d0.ToString(format);
			}
			
			public System.Random rand {
				get {
					return new System.Random(this.seed);
				}
			}
		}

		
		
		public BoundsdTests ()
		{
			System.Random rand = new System.Random("there's a rhinoceros in my soup".GetHashCode());
			AlterRollProfiles();
			//			rand = new System.Random();
			int length = Math.Max (expandedNumberOfTestItems, numberOfTestItems);
			testItemSets = new TestItemSet[length];
			for(int i=0; i<length; ++i){
				testItemSets[i] = GenerateTestItemSet(i, rand);
			}
			Boundsd.scratchGrounds();
		}

		void AlterRollProfiles ()
		{
			SetExtremeValues(boundsCenterRoller, useExtremeValuesInBoundsCenter);
			SetExtremeValues(boundsExtentsRoller, useExtremeValuesInBoundsExtents);
			SetExtremeValues(worldCoordinateRoller, useExtremeValuesInPoints);
			SetExtremeValues(rayOriginRoller, useExtremeValuesInRayOrigins);
			SetExtremeValues(rayDirectionRoller, useExtremeValuesInRayDirections);
			SetExtremeValues(valueRoller, useExtremeValuesInRawValues);


//			boundsExtentsRoller.positiveInfinityChance = 0f;
//			boundsExtentsRoller.nanChance = 0f;
		}

		void SetExtremeValues (NumberRollProfile roller, bool active)
		{
			//	So far I've verified that unity does indeed allow NaN and infinity values in its bounds
			//	Not that I'm suprised. But I still don't know why they're properties and not fields.
			if(active){
				roller.extremeValueChance = 0.1f;
				roller.nanChance = 1f;
				roller.negativeInfinityChance = 1f;
				roller.positiveInfinityChance = 1f;
				//	Including min and max values is meaningless in tests that compare the behavior
				//	float and double implementations, because the doubles simply won't exhibit the 
				//	same behavior since they won't overflow.
				//	It would be valid if and only if the double version was ALSO set to max or min value.
				//	...	 that's not hard to implement.

				roller.minValueChance = 0f;	
				roller.maxValueChance = 0f;
				roller.negativeMinValueChance = 0f;	
				roller.negativeMaxValueChance = 0f;
			} else {
				roller.extremeValueChance = 0f;
			}
		}
		
		TestItemSet GenerateTestItemSet(int index, System.Random rand){
			TestItemSet set = new TestItemSet();
			set.fq0 = GenerateRandomQuaternion(rand);
			set.fq1 = GenerateRandomQuaternion(rand);
			set.fv0 = GenerateRandomVector3(rand);
			set.fv1 = GenerateRandomVector3(rand);
			set.f0 = (float)rand.NextDouble();
			set.fb0 = GenerateRandomBounds(rand);
			set.fb1 = GenerateRandomBounds(rand);
			set.fray = GenerateRandomRay(rand);

			set.dq0 = new Quaterniond(set.fq0);
			set.dq1 = new Quaterniond(set.fq1);
			set.dv0 = new Vector3d(set.fv0);
			set.dv1 = new Vector3d(set.fv1);
			set.d0 = set.f0;
			set.db0 = new Boundsd(set.fb0);
			set.db1 = new Boundsd(set.fb1);
			set.dray = new Rayd(set.fray);
			
			set.seed = rand.Next();
			//			TestNormalizationOfSet(set);

			MatchExtremeValues(ref set);

			if(prenormalizeExtents){
				NormalizeExtents(ref set.fb0, ref set.db0);
				NormalizeExtents(ref set.fb1, ref set.db1);
			}
			return set;
		}

		void MatchExtremeValues(ref TestItemSet set){
			MatchExtremeValues(set.fb0, ref set.db0);
			MatchExtremeValues(set.fb1, ref set.db1);
			MatchExtremeValues(set.fv0, ref set.dv0);
			MatchExtremeValues(set.fv1, ref set.dv1);
			MatchExtremeValues(set.fq0, ref set.dq0);
			MatchExtremeValues(set.fq1, ref set.dq1);
			MatchExtremeValues(set.fray, ref set.dray);
			MatchExtremeValues(set.f0, ref set.d0);
		}

		void MatchExtremeValues(Bounds fb, ref Boundsd db){
			db.center = MatchExtremeValues(fb.center, db.center);
			db.extents = MatchExtremeValues(fb.extents, db.extents);
		}
		void MatchExtremeValues(Vector3 fv, ref Vector3d dv){
			MatchExtremeValues(fv.x, ref dv.x);
			MatchExtremeValues(fv.y, ref dv.y);
			MatchExtremeValues(fv.z, ref dv.z);
		}
		void MatchExtremeValues(Quaternion fv, ref Quaterniond dv){
			MatchExtremeValues(fv.x, ref dv.x);
		}
		void MatchExtremeValues(Ray fv, ref Rayd dv){
			dv.origin = MatchExtremeValues(fv.origin, dv.origin);
			dv.direction = MatchExtremeValues(fv.direction, dv.direction);
		}
		Vector3d MatchExtremeValues(Vector3 fv, Vector3d dv){
			MatchExtremeValues(fv.x, ref dv.x);
			MatchExtremeValues(fv.y, ref dv.y);
			MatchExtremeValues(fv.z, ref dv.z);
			return dv;
		}

		void MatchExtremeValues(float f, ref double d){
			if(f == float.MinValue){
				d = double.MinValue;
			} else if(f == -float.MinValue){
				d = -double.MinValue;
			} else if(f == float.MaxValue){
				d = double.MaxValue;
			} else if(f == -float.MaxValue){
				d = -double.MaxValue;
			}
		}
		double MatchExtremeValues(float f, double d){
			if(f == float.MinValue){
				return double.MinValue;
			} else if(f == -float.MinValue){
				return -double.MinValue;
			} else if(f == float.MaxValue){
				return double.MaxValue;
			} else if(f == -float.MaxValue){
				return -double.MaxValue;
			}
			return d;
		}



		
		void TestNormalizationOfSet(TestItemSet set){
			
			if(!IsNormalized(set.fq0, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.fq0 + " = " + GetSumOfSquares(set.fq0));
			}
			if(!IsNormalized(set.fq1, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.fq1 + " = " + GetSumOfSquares(set.fq1));
			}
			if(!IsNormalized(set.dq0, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.dq0 + " = " + GetSumOfSquares(set.dq0));
			}
			if(!IsNormalized(set.dq1, floatPrecision)){
				Debug.LogWarning("Not normalized: " + set.dq1 + " = " + GetSumOfSquares(set.dq1));
			}
		}
		
		float RollBoundedDouble(System.Random rand, float min, float max, float boundaryChance){
			if(rand.NextDouble() < boundaryChance){
				switch(rand.Next() % 3){
				case 0:
					return min;
				case 1:
					return max;
				default:
					return (min + max) * 0.5f;
				}
			} else {
				return (float)(rand.NextDouble() * (max - min) + min);
			}
		}

		[Test]
		public void TestTestItems (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(set.f0, set.d0);
			AssertSimilar(set.fq0, set.dq0);
			AssertSimilar(set.fq1, set.dq1);
			AssertSimilar(set.fv0, set.dv0);
			AssertSimilar(set.fv1, set.dv1);
			AssertSimilar(set.fb0, set.db0);
			AssertSimilar(set.fb1, set.db1);
		}
		double ToleranceBasisOf (Vector3d vector){
			double basis = double.MinValue;
			for(int i=0; i<3; ++i){
				double d = vector[i];
				if(d != 0d){
					if(Mathd.Abs(d) > basis){
						basis = Mathd.Abs(d);
					}
				}
			}
			return basis; 
		}
		double ToleranceBasisOf (params Vector3d[] vectors){
			double basis = double.MinValue;
			foreach(Vector3d v3 in vectors){
				double result = ToleranceBasisOf(v3);
				if(result > basis)
					basis = result;
			}
			return basis;
		}
		double ToleranceBasisOf (params double[] values){
			double basis = double.MinValue;
			foreach(double d in values){
				if(d != 0d){
					if(Mathd.Abs(d) > basis){
						basis = Mathd.Abs(d);
					}
				}
			}
			return basis;
		}
		double ToleranceBasisOf (Boundsd b)
		{
			return Mathd.Max (
				b.min.magnitude,
				b.max.magnitude,
				b.center.magnitude,
				b.size.magnitude
				);
		}

		[Test]
		public void TestMin (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(set.fb0.min, set.db0.min, ToleranceBasisOf(set.db0), "Set 0");
			AssertSimilar(set.fb1.min, set.db1.min, ToleranceBasisOf(set.db1), "Set 1");
		}
		[Test]
		public void TestMax (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(set.fb0.max, set.db0.max, ToleranceBasisOf(set.db0), "Set 0");
			AssertSimilar(set.fb1.max, set.db1.max, ToleranceBasisOf(set.db0), "Set 1");
		}
		[Test]
		public void TestCenter (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(set.fb0.center, set.db0.center, "Set 0");
			AssertSimilar(set.fb1.center, set.db1.center, "Set 1");
		}
		[Test]
		public void TestExtents (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(set.fb0.extents, set.db0.extents, "Set 0");
			AssertSimilar(set.fb1.extents, set.db1.extents, "Set 1");
		}
		[Test]
		public void TestSize (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			AssertSimilar(set.fb0.size, set.db0.size, "Set 0");
			AssertSimilar(set.fb1.size, set.db1.size, "Set 1");
		}

		[Test]
		public void TestContains (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			//	Unfortunately, sometimes a float-based bounds doesn't contain it's own min or max.
			//	This is an unfortunate side-effect of defining the bounds by its center and extents.
			//	HOWEVER. Given that most bounding boxes are going to be transformed heavily, this 
			//	representation almost certainly retains more accuraccy over multiple transformations.

			//	This actually tells me I'm on the right track with my implementation. If it was
			//	determining containment by its calculated min&max, then this would pass (or fail) 
			//	every time, not pass 10 out of 12 times. So it must be calculating the difference
			//	and comparing it to the extents.
			System.Random rand = set.rand;
			AssertContains (set.fb0.Contains(set.fb0.center), set.db0, set.db0.center, "set.db0.center");
//			Assert.AreEqual (set.fb0.Contains(set.fb0.min), set.db0.Contains(set.db0.min), "Contains(min)");
//			Assert.AreEqual (set.fb0.Contains(set.fb0.max), set.db0.Contains(set.db0.max), "Contains(max)");
//			Assert.True (set.fb0.Contains(set.fb0.center));
//			Assert.True (set.fb0.Contains(set.fb0.min));
//			Assert.True (set.fb0.Contains(set.fb0.max));
//
//			Assert.True (set.db0.Contains(set.db0.center));
//			AssertContains(set.db0, set.db0.min, "set.db0.min");
//			AssertContains(set.db0, set.db0.max, "set.db0.max");

			AssertContains(set.fb0.Contains(set.fv0), set.db0, set.dv0, "set.dv0");
			AssertContains(set.fb0.Contains(set.fv1), set.db0, set.dv1, "set.dv1");
//			Assert.AreEqual (set.fb0.Contains(set.fv0), set.db0.Contains(set.dv0), "Contains(v0=" + set.dv0 + ")");
//			Assert.AreEqual (set.fb0.Contains(set.fv1), set.db0.Contains(set.dv1), "Contains(v1=" + set.dv1 + ")");

			Bounds fb = set.fb0;
			Boundsd db = set.db0;

			Vector3[] pointsIn = RollPointsWithinBounds(rand, set.fb0, 12);
			Vector3[] pointsOut = RollPointsOutsideBounds(rand, set.fb0, 12);
			foreach(Vector3 point in pointsIn){
				AssertContains(fb.Contains(point), db, point, "");
			}
			foreach(Vector3 point in pointsOut){
				AssertContains(fb.Contains(point), db, point, "");
			}

			NormalizeExtents(ref fb, ref db);

			foreach(Vector3 point in pointsIn){
				AssertContains(fb, point, "[ALERT]");
				AssertContains(db, point, "");
			}
			foreach(Vector3 point in pointsOut){
				AssertDoesNotContain(fb, point, "[ALERT]set.fv0");
				AssertDoesNotContain(db, point, "set.dv0");
			}
		}

		void NormalizeExtents (ref Bounds fb, ref Boundsd db)
		{
			fb.extents = new Vector3(Mathf.Abs(fb.extents.x),Mathf.Abs(fb.extents.y),Mathf.Abs(fb.extents.z));
			db.extents = new Vector3d(Mathd.Abs(db.extents.x),Mathd.Abs(db.extents.y),Mathd.Abs(db.extents.z));
		}

		public static void AssertContains (Bounds bounds, Vector3 point, string info)
		{
			AssertContains (true, bounds, point, info);
		}
		public static void AssertDoesNotContain (Bounds bounds, Vector3 point, string info)
		{
			AssertContains (false, bounds, point, info);
		}
		public static void AssertContains (bool expectedResult, Bounds bounds, Vector3 point, string info)
		{
			if(bounds.Contains(point) != expectedResult){
				Assert.Fail("Bounds " + (expectedResult ? "doesn't" : "should not") + " contain point" 
				            + "\nBounds: " + StrMultiline (bounds, 2) 
				            + "\nPoint: " + Str (point) 
				            + "\nAdditional Info: " + info);
			}
		}
		public static void AssertContains (Boundsd bounds, Vector3d point, string info)
		{
			AssertContains (true, bounds, point, info);
		}
		public static void AssertDoesNotContain (Boundsd bounds, Vector3d point, string info)
		{
			AssertContains (false, bounds, point, info);
		}
		public static void AssertContains (bool expectedResult, Boundsd bounds, Vector3d point, string info)
		{
			if(bounds.Contains(point) != expectedResult){
				Assert.Fail("Bounds " + (expectedResult ? "doesn't" : "should not") + " contain point" 
				            + "\nBounds: " + StrMultiline (bounds, 2) 
				            + "\nPoint: " + Str (point) 
				            + "\nAdditional Info: " + info);
			}
		}

		[Test]
		public void TestEncapsulatePoint (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;
			
			string inputs = fb + ".Encapsulate(" + set.fv0 + ")\n"
				+ db + ".Encapsulate(" + set.dv0 + ")";

			fb.Encapsulate(set.fv0);
			db.Encapsulate(set.dv0);

			AssertSimilar(fb, db, Mathd.Max (ToleranceBasisOf(db), ToleranceBasisOf(set.db0)), inputs);
			if(ContainsPositiveInfinity(set.fv0)){
				Debug.Log ("TestEncapsulatePoint(" + testIndex + "): " + inputs + "\n= " + fb + "\n" + "= " + db);
			}
		}

		[Test]
		public void TestEncapsulateBounds (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;
			
			string inputs = fb + ".Encapsulate(" + set.fb1 + ")\n"
				+ db + ".Encapsulate(" + set.db1 + ")";

			fb.Encapsulate(set.fb1);
			db.Encapsulate(set.db1);

			AssertSimilar(fb, db, Mathd.Max (ToleranceBasisOf(db), ToleranceBasisOf(set.db0)), inputs);
//			Debug.Log ("TestEncapsulateBounds(" + testIndex + "):" + inputs + "\n= " + fb + "\n" + "= " + db);
		}

		[Test]
		public void TestExpandValue (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;
			
			string inputs = fb + ".Expand(" + set.f0 + ")\n"
				+ db + ".Expand(" + set.d0 + ")";

			fb.Expand(set.f0);
			db.Expand(set.d0);


			AssertSimilar(fb, db, Mathd.Max (ToleranceBasisOf(db), ToleranceBasisOf(set.db0)), inputs);
		}

		[Test]
		public void TestExpandVector (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;
			
			string inputs = fb + ".Expand(" + set.fv0 + ")\n"
				+ db + ".Expand(" + set.dv0 + ")";

			fb.Expand(set.fv0);
			db.Expand(set.dv0);
			
			AssertSimilar(fb, db, Mathd.Max (ToleranceBasisOf(db), ToleranceBasisOf(set.db0)), inputs);
		}

		[Test]
		public void TestIntersectRayOneArg (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;

			NormalizeExtents(ref fb, ref db);

			TestIntersectRayOneArg(fb, db, set.fray, set.dray, 0);
			TestIntersectRayOneArg(fb, db, new Ray(fb.min, fb.max - fb.min), new Rayd(db.min, db.max - db.min), 1);
			TestIntersectRayOneArg(fb, db, new Ray(fb.max + fb.size, fb.max - fb.min), new Rayd(db.max + db.size, db.max - db.min), -1);

		}

		void TestIntersectRayOneArg (Bounds fb, Boundsd db, Ray fray, Rayd dray, int knownResult)
		{
			bool fbool = fb.IntersectRay(fray);
			bool dbool = db.IntersectRay(dray);

			Assert.AreEqual(fbool, dbool);
			if(knownResult == 1){
				Assert.True(fbool);
				Assert.True(dbool);
			} else if(knownResult == -1){
				Assert.False(fbool);
				Assert.False(dbool);
			}
		}

		[Test]
		public void TestIntersectRayTwoArg (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;

//			TestIntersectRayTwoArg(fb, db, set.fray, 0);
//			TestIntersectRayTwoArg(fb, db, new Ray(fb.min + fb.size * 0.1f, fb.max - fb.min), 0);
//			TestIntersectRayTwoArg(fb, db, new Ray(fb.max + fb.size, fb.max - fb.min), 0);

			NormalizeExtents(ref fb, ref db);

//			TestIntersectRayTwoArg(fb, db, set.fray, set.dray, 0);
//			TestIntersectRayTwoArg(fb, db, new Ray(fb.min + fb.size * 0.1f, fb.max - fb.min), new Rayd(db.min + db.size * 0.1d, db.max - db.min), 1);
//			TestIntersectRayTwoArg(fb, db, new Ray(fb.max + fb.size, fb.max - fb.min), new Rayd(db.max + db.size, db.max - db.min), -1);
			
			TestIntersectRayTwoArg(fb, db, set.fray, 0);
			TestIntersectRayTwoArg(fb, db, new Ray(fb.min + fb.size * 0.1f, fb.max - fb.min), 1);
			TestIntersectRayTwoArg(fb, db, new Ray(fb.max + fb.size, fb.max - fb.min), -1);
		}
		void TestIntersectRayTwoArg (Bounds fb, Boundsd db, Ray fray, int knownResult){
			TestIntersectRayTwoArg(fb, db, fray, new Rayd(fray), knownResult);
		}
		void TestIntersectRayTwoArg (Bounds fb, Boundsd db, Ray fray, Rayd dray, int knownResult)
		{
			string inputs = "Bounds: " + StrMultiline(db,2) + "\nRay: " + Str (dray) + "\nknownResult=" + knownResult;
			float f;
			double d;
			bool fbool = fb.IntersectRay(fray, out f);
			bool dbool = db.IntersectRay(dray, out d);

			Vector3 fpoint = fray.GetPoint(f);
			Vector3d dpoint = dray.GetPoint(d);

			inputs += "\nDistances: " + f + ", " + d + 
				"\nPoints: " + fpoint.ToString("G5") + ", " + dpoint.ToString("G5") +
				"\nBools: " + fbool + ", " + dbool;

			//			Debug.Log ("Tolerance Basic(" + f + "," + d +")= " + toleranceBasis);
			AssertSimilar(fpoint, dpoint, 4*ToleranceBasisOf(fpoint, dpoint), inputs + "\npoints");
			AssertSimilar(f, d, inputs + "\ndistance", 4d * ToleranceBasisOf(f, d));

			if(knownResult == 1){
				Assert.True(fbool,"[ALERT:]\n" + inputs + "\nFailed to intersect!");
				Assert.True(dbool, inputs + "\nFailed to intersect!");
			} else if(knownResult == -1){
				Assert.False(fbool,"[ALERT:]\n" + inputs + "\nShould not intersect!");
				Assert.False(dbool, inputs + "\nShould not intersect!");
			}
			Assert.AreEqual(fbool, dbool, inputs);
		}

		int intersections = 0;
		[Test]
		public void TestIntersects (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];

			intersections += IntersectionTest(set.fb0, set.db0, set.fb1, set.db1);
		}


		int IntersectionTest (Bounds fb0, Boundsd db0, Bounds fb1, Boundsd db1)
		{
			WeakIntersectionTest(fb0,db0,fb1,db1);
			
			//	After normalization
			NormalizeExtents(ref fb0, ref db0);
			NormalizeExtents(ref fb1, ref db1);

			if(ContainsNaN(fb0,fb1)){
				WeakIntersectionTest(fb0,db0,fb1,db1);
			} else {
				StrongIntersectionTest(fb0,db0,fb1,db1);
			}
			return db0.Intersects(db1) ? 1 : 0;
		}
		void WeakIntersectionTest(Bounds fb0, Boundsd db0, Bounds fb1, Boundsd db1){
			string input0 =  StrMultiline(db0, 2);
			string input1 =  StrMultiline(db1, 2);
			string inputs = "Bounds 0: " + input0 + "\nBounds 1: " + input1;
			Assert.AreEqual(fb0.Intersects(fb0), db0.Intersects(db0), "Bounds: " + input0);
			Assert.AreEqual(fb1.Intersects(fb1), db1.Intersects(db1), "Bounds: " + input1);
			Assert.AreEqual(fb0.Intersects(fb1), db0.Intersects(db1), inputs);
			Assert.AreEqual(fb1.Intersects(fb0), db1.Intersects(db0), "[ALERT:INVERSE OF PASSED OPERATION]\n" + inputs);
		}
		void StrongIntersectionTest(Bounds fb0, Boundsd db0, Bounds fb1, Boundsd db1){
			string input0 =  StrMultiline(db0, 2);
			string input1 =  StrMultiline(db1, 2);
			string inputs = "Bounds 0: " + input0 + "\nBounds 1: " + input1;

			Assert.True (fb0.Intersects(fb0), "[ALERT]Bounds: " + input0);
			Assert.True (db0.Intersects(db0), "Bounds: " + input0);
			Assert.True (fb1.Intersects(fb1), "[ALERT]Bounds: " + input1);
			Assert.True (db1.Intersects(db1), "Bounds: " + input1);
			
			Assert.AreEqual(fb0.Intersects(fb1), db0.Intersects(db1), inputs);
			Assert.AreEqual(fb1.Intersects(fb0), db1.Intersects(db0), "[ALERT:INVERSE OF PASSED OPERATION]\n" + inputs);
		}

		[Test]
		public void TestIntersectsOverAllTestSets (){
			int numTests = 0, numPassed = 0;
			for(int i=0; i<numberOfTestItems-1; ++i){
				TestItemSet set0 = testItemSets[i];
				for(int j=i+1; j<numberOfTestItems; ++j){
					TestItemSet set1 = testItemSets[j];
					numPassed += IntersectionTest(set0.fb0, set0.db0, set1.fb0, set1.db0);
					numPassed += IntersectionTest(set0.fb0, set0.db0, set1.fb1, set1.db1);
					numPassed += IntersectionTest(set0.fb1, set0.db1, set1.fb1, set1.db1);
					numPassed += IntersectionTest(set0.fb1, set0.db1, set1.fb0, set1.db0);
					numTests += 4;
				}
			}
			//	613/7080 using expanded number of test items. All tests passed. Pretty solid implementation, then.
//			Debug.Log ("Tests Intersected: " + numPassed  + "/" + numTests + "\nNormal: " + intersections + "/" + numberOfTestItems);
		}

		[Test]
		public void TestSetMinMax (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;

			fb.SetMinMax(set.fv0, set.fv1);
			db.SetMinMax(set.dv0, set.dv1);
			
			AssertSimilar (fb, db, ToleranceBasisOf(db));
		}

		[Test]
		public void TestSqrDistance (
			[NUnit.Framework.Range (0,numberOfTestItems-1)] int testIndex
			){
			TestItemSet set = testItemSets[testIndex];
			
			Bounds fb = set.fb0;
			Boundsd db = set.db0;

			Vector3d diff, nonNaiveDiff;

			float f = fb.SqrDistance(set.fv0);
			double d = db.SqrDistance(set.dv0, out diff);
			double nonNaive = db.SqrDistanceNotNaively(set.dv0, out nonNaiveDiff);


			string inputs = "Bounds: " + StrMultiline(set.db0, 2) + "\nPoint: " + Str (set.dv0);

			inputs += "\nDiff: " + diff + " = " + diff.sqrMagnitude + " => " + diff.magnitude;
			inputs += "\nF Result: " + f + " => " + Mathf.Sqrt(f);
			double tolerance = 2 * Math.Min (Math.Abs(f), Math.Abs(d)) ;
			if(Math.Abs (f - d) >= tolerance * floatPrecision){
				inputs += seekOutDiffMatch(f, diff, nonNaiveDiff);
			}
			AssertSimilar (f, d, inputs + "\nsqrDistance", tolerance);
			Debug.Log ("TestSqrDistance(" + testIndex + "): \n" + inputs);
		}

		string seekOutDiffMatch (float f, Vector3d diff, Vector3d nonNaiveDiff)
		{
			StringBuilder b = new StringBuilder();
//			AppendDiffMatch(b, f, diff.x, diff.y, diff.z);
//			AppendDiffMatch(b, f, 0, diff.y, diff.z);
//			AppendDiffMatch(b, f, diff.x, 0, diff.z);
//			AppendDiffMatch(b, f, diff.x, diff.y, 0);
//			AppendDiffMatch(b, f, nonNaiveDiff.x, diff.y, diff.z);
//			AppendDiffMatch(b, f, diff.x, nonNaiveDiff.y, diff.z);
//			AppendDiffMatch(b, f, diff.x, diff.y, nonNaiveDiff.z);

			
			AppendDiffMatch(b, f, Mathd.Sqrt(f - diff.y * diff.y - diff.z * diff.z), diff.y, diff.z);
			AppendDiffMatch(b, f, diff.x, Mathd.Sqrt(f - diff.x * diff.x - diff.z * diff.z), diff.z);
			AppendDiffMatch(b, f, diff.x, diff.y, Mathd.Sqrt(f - diff.y * diff.y - diff.x * diff.x));
//			AppendDiffMatch(b, f, diff.x, diff.y, diff.z);

			return b.ToString ();
		}

		void AppendDiffMatch (StringBuilder b, float f, double x, double y, double z)
		{
			Vector3d diff = new Vector3d(x,y,z);
			double d = diff.sqrMagnitude;
			double tolerance = 2 * Math.Min (Math.Abs(f), Math.Abs(d)) * floatPrecision;

			b.Append("\n").Append (diff).Append(" => ").Append(d).Append(" => ").Append(Mathd.Abs(f - d) / tolerance);
		}

	}
}

