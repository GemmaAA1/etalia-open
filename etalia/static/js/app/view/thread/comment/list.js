define([
    'app',
    'text!app/templates/thread/comment/list.hbs',
    'app/view/thread/comment/thumb',
    'app/view/thread/comment/form'
], function (App, template) {

    App.View.Thread.CommentList = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-comments',

        template: App.Handlebars.compile(template),

        post: null,

        events: {
            "click .tread-comment-add": "onAddCommentClick"
        },

        initialize: function (options) {
            if (!options.post) {
                throw 'Expected instance of App.Model.Post';
            }

            this.post = options.post;

            this.listenTo(this.model, "sync", this.render);
        },

        onAddCommentClick: function(e) {
            e.preventDefault();

            var that = this,
                $link = this.$('.tread-comment-add').hide();

            var form = App.View.Thread.CommentForm.create({},{
                $target: this.$('[data-comment-add-form]')
            });

            var cancel = function() {
                form.$el.replaceWith('<div data-comment-add-form></div>');
                form.remove();
                $link.show();
            };

            this.listenToOnce(form, 'cancel', cancel);
            this.listenToOnce(form, 'validation_success', function() {
                form.model.save({
                    user: App.getCurrentUser(),
                    post: that.post
                },{
                    success: cancel
                });
            });
        },

        render: function() {
            App.log('CommentListView::render');

            this.$el.html(this.template({
                is_member: this.post.get('thread').isMember(App.getCurrentUser())
            }));

            var that = this,
                $list = this.$('.thread-comments-list');
            this.model.each(function(comment) {
                that.pushSubView(
                    App.View.Thread.CommentThumb.create({
                        model: comment
                    }, {
                        $target: $list,
                        append: true
                    })
                );
            });

            return this;
        }
    });

    App.View.Thread.CommentList.create = function(options, createOptions) {
        options = options || {};
        if (!options.post) {
            throw 'options.post is expected to be an instace of App.Model.Post';
        }
        options.model = options.post.get('comments');

        var list = new App.View.Thread.CommentList(options);

        if (createOptions) {
            App.View.create(list, createOptions);
        }
        return list;
    };

    return App.View.Thread.CommentList;
});
