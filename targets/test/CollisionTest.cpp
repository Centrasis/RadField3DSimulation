#include "Collisions.h"
#include <iostream>
#include "gtest/gtest.h"

namespace {
	Collisions::Box box(glm::vec3(-1.f), glm::vec3(1.f));
	Collisions::Triangle triangle(glm::vec3(-1.f), glm::vec3(1.f, -1.f, -1.f), glm::vec3(1.f, 1.f, -1.f));

	TEST(TriangleDissectionTest, CenteredLine) {
		EXPECT_TRUE(triangle.intersect(Collisions::Line(glm::vec3(0.5f, -0.5f, 2.f), glm::vec3(0.5f, -0.5f, -2.f))));
	};

	TEST(TriangleDissectionTest, CenteredLineMirrored) {
		EXPECT_TRUE(triangle.intersect(Collisions::Line(glm::vec3(0.5f, -0.5f, -2.f), glm::vec3(0.5f, -0.5f, 2.f))));
	};

	TEST(TriangleDissectionTest, CenteredPoint) {
		EXPECT_TRUE(triangle.intersect(Collisions::Line(glm::vec3(0.5f, -0.5f, -1.f), glm::vec3(0.5f, -0.5f, -1.f))));
	};

	TEST(TriangleMissTest, AlignedLine) {
		EXPECT_FALSE(triangle.intersect(Collisions::Line(glm::vec3(-2.f, 2.f, 0.f), glm::vec3(2.f, 2.f, 0.f))));
	};

	TEST(TriangleMissTest, ShortLine) {
		EXPECT_FALSE(triangle.intersect(Collisions::Line(glm::vec3(0.5f, -0.5f, 0.f), glm::vec3(0.5f, -0.5f, -0.9f))));
	};

	TEST(TriangleMissTest, ShortLineMirrored) {
		EXPECT_FALSE(triangle.intersect(Collisions::Line(glm::vec3(0.5f, -0.5f, -2.f), glm::vec3(0.5f, -0.5f, -1.01f))));
	};

	TEST(TriangleMissTest, Point) {
		EXPECT_FALSE(triangle.intersect(Collisions::Line(glm::vec3(0.f, 2.f, 0.f), glm::vec3(0.f, 2.f, 0.f))));
	};

	TEST(BoxDissectionTest, CenteredLine) {
		EXPECT_TRUE(box.intersect(Collisions::Line(glm::vec3(-2.f, 0.f, 0.f), glm::vec3(2.f, 0.f, 0.f))));
	};

	TEST(BoxDissectionTest, CrossingLine) {
		EXPECT_TRUE(box.intersect(Collisions::Line(glm::vec3(-2.f, 0.f, 0.f), glm::vec3(2.f, 1.f, 1.f))));
	};

	TEST(BoxDissectionTest, CenteredPoint) {
		EXPECT_TRUE(box.intersect(Collisions::Line(glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 0.f))));
	};

	TEST(BoxMissTest, AlignedLine) {
		EXPECT_FALSE(box.intersect(Collisions::Line(glm::vec3(-2.f, 2.f, 0.f), glm::vec3(2.f, 2.f, 0.f))));
	};

	TEST(BoxMissTest, Point) {
		EXPECT_FALSE(box.intersect(Collisions::Line(glm::vec3(0.f, 2.f, 0.f), glm::vec3(0.f, 2.f, 0.f))));
	};
};
