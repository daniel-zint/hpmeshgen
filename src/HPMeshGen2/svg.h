#pragma once

#include <fstream>
#include <experimental/filesystem>
#include "MeshHeader.h"

namespace svg{
	inline void svgRect(std::ofstream& ofs, float x1, float y1, float x2, float y2, std::string fill = "white", std::string stroke = "black", float stroke_width = 1) {
		ofs << "<rect x=\"" << x1 << "\" y=\"" << y1 << "\" width=\"" << x2 - x1 << "\" height=\"" << y2 - y1 << "\" fill=\"" << fill << "\" stroke=\"" << stroke << "\" stroke-width=\"" << stroke_width << "px\"/>" << std::endl;
	}

	inline void svgLine(std::ofstream& ofs, float x1, float y1, float x2, float y2, std::string stroke = "black", float stroke_width = 1) {
		ofs << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" stroke=\"" << stroke << "\" stroke-width=\"" << stroke_width << "px\" stroke-linecap=\"round\"/>" << std::endl;
	}

	inline void svgLine(std::ofstream& ofs, const TriMesh::Point& p1, const TriMesh::Point& p2, std::string stroke = "black", float stroke_width = 1) {
		svgLine(ofs, p1[0], p1[1], p2[0], p2[1], stroke, stroke_width);
	}
	
	inline void svgCircle(std::ofstream& ofs, float x, float y, float r = 1, std::string fill = "white", std::string stroke = "black", float stroke_width = 1) {
		ofs << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << r << "\" fill=\"" << fill << "\" stroke=\"" << stroke << "\" stroke-width=\"" << stroke_width << "px\"/>" << std::endl;
	}

	inline void svgCircle(std::ofstream& ofs, TriMesh::Point p, float r = 1, std::string fill = "white", std::string stroke = "black", float stroke_width = 1) {
		svgCircle(ofs, p[0], p[1], r, fill, stroke, stroke_width);
	}

	inline void svgPolygon(std::ofstream& ofs, const std::vector<float>& x, const std::vector<float>& y, const std::string& style = "fill:rgba(255,0,0,0.5);stroke:black;stroke-width:0") {
		assert(x.size() == y.size());

		ofs << "<polygon points = \"";

		for (size_t i = 0; i < x.size(); ++i) {
			ofs << x[i] << "," << y[i] << " ";
		}
		ofs << "\" style = \"" << style << "\" />" << std::endl;
	}

	inline void svgPolygon(std::ofstream& ofs, const std::vector<TriMesh::Point>& p, const std::string& style = "fill:rgba(255,0,0,0.5);stroke:black;stroke-width:0") {
		std::vector<float> x(p.size());
		std::vector<float> y(p.size());

		for (size_t i = 0; i < p.size(); ++i) {
			x[i] = p[i][0];
			y[i] = p[i][1];
		}

		svgPolygon(ofs, x, y, style);
	}



	template <typename T>
	inline void printInteriorEdges(std::ofstream& ofs, T& mesh) {
		for (auto eh_it = mesh.edges_begin(); eh_it != mesh.edges_end(); ++eh_it) {
			if (mesh.is_boundary(*eh_it))
				continue;

			typename T::HalfedgeHandle heh = mesh.halfedge_handle(*eh_it, 0);

			typename T::Point p1 = mesh.point(mesh.from_vertex_handle(heh));
			typename T::Point p2 = mesh.point(mesh.to_vertex_handle(heh));

			svgLine(ofs, p1, p2);
		}
	}

	template <typename T>
	inline void printBoundaryEdges(std::ofstream& ofs, T& mesh) {
		for (auto eh_it = mesh.edges_begin(); eh_it != mesh.edges_end(); ++eh_it) {
			if (!mesh.is_boundary(*eh_it))
				continue;

			typename T::HalfedgeHandle heh = mesh.halfedge_handle(*eh_it, 0);

			typename T::Point p1 = mesh.point(mesh.from_vertex_handle(heh));
			typename T::Point p2 = mesh.point(mesh.to_vertex_handle(heh));

			svgLine(ofs, p1, p2, "red", 2);
		}
	}

	template <typename T> inline void printVertices(std::ofstream& ofs, T& mesh) {
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
			typename T::Point p = mesh.point(*v_it);

			svgCircle(ofs, p, 1, "black", "black", 1);
		}
	}

	template <typename T> inline void printFace(std::ofstream& ofs, T& mesh, uint fID, const std::string& style = "fill:rgba(255,0,0,0.5);stroke:black;stroke-width:0") {
		typename T::FaceHandle fh = mesh.face_handle(fID);

		std::vector<typename T::Point>p;

		for (auto fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
			p.push_back(mesh.point(*fv_it));
		}

		svgPolygon(ofs, p, style);
	}

	template <typename T> void printMeshSVG(const std::experimental::filesystem::path& filename, T& meshInput, const bool printValence = false) {
		std::ofstream ofs(filename);

		T mesh = meshInput;

		// get min/max
		float x_min = std::numeric_limits<float>::max();
		float y_min = std::numeric_limits<float>::max();
		float x_max = -std::numeric_limits<float>::max();
		float y_max = -std::numeric_limits<float>::max();

		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
			typename T::Point p = mesh.point(*v_it);

			x_min = std::min(x_min, p[0]);
			x_max = std::max(x_max, p[0]);
			y_min = std::min(y_min, p[1]);
			y_max = std::max(y_max, p[1]);
		}

		float x_dim = x_max - x_min;
		float y_dim = y_max - y_min;

		int xpx = 1000;
		int ypx = (int)(xpx * y_dim / x_dim);

		// norm all vertices
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
			typename T::Point p = mesh.point(*v_it);
			p[0] = xpx * (p[0] - x_min) / x_dim;
			p[1] = ypx - ypx * (p[1] - y_min) / y_dim;
			p[2] = 0;

			mesh.set_point(*v_it, p);
		}

		// header
		ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
		ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;
		ofs << "version=\"1.1\" baseProfile=\"full\"" << std::endl;
		ofs << "width=\"" << xpx << "px\" height=\"" << ypx << "px\" viewBox=\"0 0 " << xpx << " " << ypx << "\">" << std::endl;
		ofs << "" << std::endl;

		printInteriorEdges(ofs, mesh);
		printBoundaryEdges(ofs, mesh);

		if (printValence) {
			for (auto vh_it = mesh.vertices_begin(); vh_it != mesh.vertices_end(); ++vh_it) {
				TriMesh::Point p = mesh.point(*vh_it);

				int valence;
				if (mesh.is_boundary(*vh_it))
					valence = 4;
				else
					valence = 6;

				int n_neighbors = 0;
				for (auto vv_it = mesh.vv_iter(*vh_it); vv_it.is_valid(); ++vv_it) ++n_neighbors;

				if (n_neighbors > valence)
					svg::svgCircle(ofs, p, 10.f, "blue", "black", 1.f);
				if (n_neighbors < valence)
					svg::svgCircle(ofs, p, 10.f, "orange", "black", 1.f);
			}
		}

		ofs << "</svg>" << std::endl;

		ofs.close();
	}
}