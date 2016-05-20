define([
    'app',
    'text!app/templates/ui/filters.hbs',
    'text!app/templates/ui/filters-group.hbs',
    'text!app/templates/ui/filters-entry.hbs',
    'app/model/ui/filters'
], function (App, template, groupTemplate, entryTemplate) {

    App.View.Ui = App.View.Ui || {};

    App.View.Ui.FilterEntry = App.Backbone.View.extend({
        tagName: 'li',

        template: App.Handlebars.compile(entryTemplate),

        groupView: null,

        events: {
            "click a": "onToggleClick"
        },

        initialize: function (options) {
            this.groupView = options.groupView;

            this.listenTo(this.model, 'change', this.render);
        },

        onToggleClick: function(e) {
            e.preventDefault();

            this.model.set('active', !this.model.get('active'));

            this.groupView.filtersView.trigger('context-change');
        },

        render: function () {
            this.$el.html(this.template(this.model.attributes));

            return this;
        }
    });

    App.View.Ui.FilterGroup = App.Backbone.View.extend({
        tagName: 'div',
        className: 'filter-group',

        template: App.Handlebars.compile(groupTemplate),

        filtersView: null,
        shownCount: 10,

        events: {
            "click a.filter-toggle": "onToggleClick",
            "click a.filter-more": "onMoreClick"
        },

        initialize: function (options) {
            this.filtersView = options.filtersView;
        },

        onToggleClick: function(e) {
            e.preventDefault();

            if (this.$el.hasClass('active')) {
                this.reduce();
            } else {
                this.expand();
            }
        },

        onMoreClick: function(e) {
            e.preventDefault();

            this.shownCount += 10;

            this.applyFiltersVisibility();
        },

        expand: function() {
            this.$el.addClass('active');
            this.$('.collapse').addClass('in');

            return this;
        },

        reduce: function() {
            this.$el.removeClass('active');
            this.$('.collapse').removeClass('in');

            return this;
        },

        applyFiltersVisibility: function() {
            var $filters = this.$('li');

            if (this.shownCount >= $filters.length) {
                $filters.show();
                this.$('.filter-more').hide();
            } else {
                $filters.filter(':lt(' + this.shownCount + ')').show();
                $filters.filter(':gt(' + (this.shownCount-1) + ')').hide();
                this.$('.filter-more').show();
            }

            return this;
        },

        render: function () {
            this.$el.html(this.template(this.model.attributes));

            var that = this,
                $list = this.$('ul');

            this.model.get('entries').each(function(entry) {
                var entryView = new App.View.Ui.FilterEntry({
                    groupView: that,
                    model: entry
                });

                $list.append(entryView.render().$el);
                that.pushSubView(entryView);
            });

            this.applyFiltersVisibility();

            return this;
        }
    });

    App.View.Ui.Filters = App.Backbone.View.extend({
        tagName: 'div',

        attributes: {
            id: 'filter-flap',
            'class': 'flap right'
        },

        template: App.Handlebars.compile(template),

        setGroups: function(groups) {
            this.clear();
            this.model = groups;
            this.render();

            return this;
        },

        getContext: function() {
            var result = {};

            this.model.each(function (group) {
                result[group.get('name')] = group.get('entries')
                    .filter(function (entry) {
                        return entry.get('active');
                    })
                    .map(function(entry) {
                        return entry.get('id');
                    });
            });

            return result;
        },

        clear: function() {
            if (this.subViews) {
                App._.each(this.subViews, function (subView) {
                    subView.remove();
                });
                this.subViews = [];
            }

            return this;
        },

        load: function(url, data) {
            var that = this,
                groups = new App.Model.Ui.FilterGroups();

            groups.fetch({
                    url: url,
                    data: data
                })
                .done(function() {
                    that.setGroups(groups);
                    that.trigger('loaded');
                })
                .fail(function() {
                    App.log(arguments);
                    throw 'Failed to fetch "' + url + '" filters.';
                });
        },

        render: function () {
            this.$el.html(this.template());

            var that = this,
                $list = this.$('.flap-inner');

            this.model.each(function(group, i) {
                var groupView = new App.View.Ui.FilterGroup({
                    filtersView: that,
                    model: group
                });

                $list.append(groupView.render().$el);

                if (i == 0) {
                    groupView.expand();
                }

                that.pushSubView(groupView);
            });

            return this;
        }
    });

    App.View.Ui.Filters.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = new App.Model.Ui.FilterGroups();
        }

        var view = new App.View.Ui.Filters(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Ui.Filters;
});
