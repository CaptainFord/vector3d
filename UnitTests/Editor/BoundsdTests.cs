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
		const int numberOfTestItems = 12;
		const int expandedNumberOfTestItems = 60;
		
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
			//			rand = new System.Random();
			testItemSets = new TestItemSet[expandedNumberOfTestItems];
			for(int i=0; i<expandedNumberOfTestItems; ++i){
				testItemSets[i] = GenerateTestItemSet(i, rand);
			}
			Boundsd.scratchGrounds();
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
			set.fray = new Ray(GenerateRandomVector3(rand), GenerateRandomVector3(rand));

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
			return set;
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
			Assert.AreEqual (set.fb0.Contains(set.fb0.center), set.db0.Contains(set.db0.center), "Contains(center)");
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

			NormalizeExtents(ref fb, ref db);

			TestIntersectRayTwoArg(fb, db, set.fray, set.dray, 0);
			TestIntersectRayTwoArg(fb, db, new Ray(fb.min + fb.size * 0.1f, fb.max - fb.min), new Rayd(db.min + db.size * 0.1d, db.max - db.min), 1);
			TestIntersectRayTwoArg(fb, db, new Ray(fb.max + fb.size, fb.max - fb.min), new Rayd(db.max + db.size, db.max - db.min), -1);
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
			AssertSimilar(f, d, inputs + "\ndistance", 50d * ToleranceBasisOf(f, d));
			AssertSimilar(fpoint, dpoint, 4*ToleranceBasisOf(fpoint, dpoint), inputs + "\npoints");

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
			string input0 =  StrMultiline(db0, 2);
			string input1 =  StrMultiline(db1, 2);
			string inputs = "Bounds 0: " + input0 + "\nBounds 1: " + input1;
			Assert.AreEqual(fb0.Intersects(fb0), db0.Intersects(db0), "Bounds: " + input0);
			Assert.AreEqual(fb1.Intersects(fb1), db1.Intersects(db1), "Bounds: " + input1);
			Assert.AreEqual(fb0.Intersects(fb1), db0.Intersects(db1), inputs);
			Assert.AreEqual(fb1.Intersects(fb0), db1.Intersects(db0), "[ALERT:INVERSE OF PASSED OPERATION]\n" + inputs);
			
			//	After normalization
			NormalizeExtents(ref fb0, ref db0);
			NormalizeExtents(ref fb1, ref db1);
			
			Assert.True (fb0.Intersects(fb0), "[ALERT]Bounds: " + input0);
			Assert.True (db0.Intersects(db0), "Bounds: " + input0);
			Assert.True (fb1.Intersects(fb1), "[ALERT]Bounds: " + input1);
			Assert.True (db1.Intersects(db1), "Bounds: " + input1);
			
			Assert.AreEqual(fb0.Intersects(fb1), db0.Intersects(db1), inputs);
			Assert.AreEqual(fb1.Intersects(fb0), db1.Intersects(db0), "[ALERT:INVERSE OF PASSED OPERATION]\n" + inputs);
			return db0.Intersects(db1) ? 1 : 0;
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

			float f = fb.SqrDistance(set.fv0);
			double d = db.SqrDistance(set.dv0);

			string inputs = "Bounds: " + StrMultiline(set.db0, 2) + "\nPoint: " + Str (set.dv0);
			
			AssertSimilar (f, d, inputs + "\nsqrDistance", Math.Min (2d*f, 2d*d));
		}
	}
}

