define([
    'app',
    'text!app/templates/thread/list.html',
    'app/view/detail',
    'app/view/thread/detail',
    'app/view/thread/thumb',
    'app/collection/thread/thread'
], function (App, template) {

    var defaults = {
        query: {
            view: 'nested'
        }
    };

    App.View.Thread = App.View.Thread || {};

    return App.View.Thread.List = App.Backbone.View.extend({

        //tagName: 'div',
        //className: 'thumb-list',
        template: App.Handlebars.compile(template),

        detailView: null,
        $thumbList: null,

        events: {
            "click #thread-create-modal": "onCreateModalClick",
            "click #thread-next-page": "onNextPageClick"
        },

        initialize: function (options) {
            App.defaults(options, defaults);

            this.collection = new App.Collection.Threads({
                query: options.query
            });

            this.listenTo(this.collection.fullCollection, "reset", this.onCollectionReset);
            this.listenTo(this.collection.fullCollection, "add", this.onCollectionAdd);
            this.listenTo(this.collection.fullCollection, "remove", this.onCollectionRemove);

            this.listenTo(this, "model:remove", this.onModelRemove);
            this.listenTo(this, "model:detail", this.openDetail);
        },

        openDetail: function (model) {
            App.log('ThreadListView::onModelDetail');

            var modelIndex = this.collection.fullCollection.indexOf(model);
            if (0 > modelIndex) {
                throw 'Unexpected index';
            }

            var detailModel = new App.Model.Detail({
                list: this,
                prev: 0 < modelIndex ? this.collection.fullCollection.at(modelIndex - 1) : null,
                next: this.collection.fullCollection.at(modelIndex + 1),
                view: new App.View.Thread.Detail({
                    model: model
                })
            });

            if (this.detailView) {
                this.detailView.model.get('view').remove();
                this.detailView.model.destroy();
                this.detailView.model = detailModel;
            } else {
                this.detailView = new App.View.Detail({
                    model: detailModel
                });
                this.listenTo(this.detailView, "detail:prev", this.openDetail);
                this.listenTo(this.detailView, "detail:next", this.openDetail);
            }

            var that = this;
            model
                .fetch({data: {view: 'nested'}})
                .done(function() {
                    that.detailView.render();
                });
        },

        onModelRemove: function (model) {
            App.log('ThreadListView::onModelRemove');

            this.collection.remove(model);
            this.collection.fullCollection.remove(model);
        },

        onCollectionAdd: function (model) {
            App.log('ThreadListView::onCollectionAdd');

            // Render the thumb
            var thumbView = new App.View.Thread.Thumb({
                    id: 'thread-thumb-' + model.get('id'),
                    model: model,
                    list: this
                });

            // Append the thumb
            thumbView
                .render()
                .$el.appendTo(this.$thumbList);
        },

        onCollectionRemove: function (model) {
            App.log('ThreadListView::onCollectionRemove');
            this.$('#thread-thumb-' + model.get('id')).remove();
        },

        onCollectionReset: function () {
            if (this.collection.fullCollection.length == 0) {
                this.render();
            }
        },

        onNextPageClick: function (e) {
            e.preventDefault();

            var that = this;
            if (this.collection.hasNextPage()) {
                this.collection.getNextPage(null).then(function() {
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                });
            }
        },

        onCreateModalClick: function(e) {
            e.preventDefault();

            var that = this,
                form = App.View.Thread.CreateForm.create(),
                modal = new App.View.Modal({
                    title: 'Start a new thread',
                    content: form,
                    footer: false
                });

            form.once('validation_success', function () {
                form.model.save(null, {
                    wait: true,
                    success: function () {
                        that.collection.add(form.model, {at: 0});
                        modal.close();
                        that.openDetail(form.model);
                    },
                    error: function () {
                        // TODO
                    }
                });
            });

            form.once('cancel', function () {
                modal.close();
            });

            modal.once('hidden', function () {
                form = null;
                modal = null;
            });

            modal.render();
        },

        render: function () {
            App.log('ThreadListView::render');

            var that = this;

            this.$el.html(this.template({}));

            this.$thumbList = this.$('.thumb-list');

            this.collection.fetch()
                .then(function() {
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                });

            return this;
        }
    });
});
